use biofile::bed::{Bed, BedDataLine, BedWriter, Chrom};
use math::{
    interval::{traits::Interval, I64Interval},
    iter::{AggregateOp, CommonRefinementZip, IntoBinnedIntervalIter},
    partition::integer_interval_map::IntegerIntervalMap,
};
use std::collections::{HashMap, HashSet};

type Coeficient = f64;
type Value = f64;

pub struct TrackZipper {
    weighted_bed_files: Vec<(Coeficient, Bed)>,
    list_of_chrom_interval_maps: Vec<HashMap<Chrom, IntegerIntervalMap<Value>>>,
}

pub struct LinearTrackMixer {
    chrom_to_binned_values:
        HashMap<Chrom, Vec<(I64Interval, Vec<Option<Value>>)>>,
    weights: Vec<Coeficient>,
}

pub struct LinearTrackMixerIter<'a> {
    binned_values_iter: std::slice::Iter<'a, (I64Interval, Vec<Option<Value>>)>,
    weights: &'a Vec<Coeficient>,
}

impl TrackZipper {
    pub fn new(
        weighted_bed_files: Vec<(Coeficient, Bed)>,
        exclude_track_filepath: Option<&str>,
    ) -> Result<Self, biofile::error::Error> {
        let exclude = if let Some(path) = exclude_track_filepath {
            Some(Bed::new(path).get_chrom_to_intervals())
        } else {
            None
        };

        let list_of_chrom_interval_maps: Vec<
            HashMap<Chrom, IntegerIntervalMap<Value>>,
        > = weighted_bed_files
            .iter()
            .map(|(_weight, bed)| {
                bed.get_chrom_to_interval_to_val(exclude.as_ref())
            })
            .collect::<Result<
                Vec<HashMap<Chrom, IntegerIntervalMap<Value>>>,
                biofile::error::Error,
            >>()?;

        Ok(TrackZipper {
            weighted_bed_files,
            list_of_chrom_interval_maps,
        })
    }

    pub fn to_linear_track_mixer(
        &self,
        target_chroms: &HashSet<Chrom>,
        bin_size: i64,
    ) -> Result<LinearTrackMixer, biofile::error::Error> {
        let chrom_to_binned_values: HashMap<
            Chrom,
            Vec<(I64Interval, Vec<Option<Value>>)>,
        > = self.chrom_to_binned_zipped_values(target_chroms, bin_size)?;

        let weights: Vec<Coeficient> =
            self.weighted_bed_files.iter().map(|(w, _)| *w).collect();

        Ok(LinearTrackMixer {
            chrom_to_binned_values,
            weights,
        })
    }

    pub fn chrom_to_binned_zipped_values(
        &self,
        target_chroms: &HashSet<Chrom>,
        bin_size: i64,
    ) -> Result<
        HashMap<Chrom, Vec<(I64Interval, Vec<Option<Value>>)>>,
        biofile::error::Error,
    > {
        let empty_interval_map = IntegerIntervalMap::new();
        let union_zipped_chrom_interval_maps: HashMap<
            Chrom,
            Vec<&IntegerIntervalMap<Value>>,
        > = crate::util::get_union_zipped_chrom_interval_maps(
            self.list_of_chrom_interval_maps.iter().collect(),
            target_chroms,
            &empty_interval_map,
        );
        let chroms =
            crate::util::get_sorted_keys(&union_zipped_chrom_interval_maps);

        Ok(chroms
            .into_iter()
            .map(|chrom| {
                let interval_maps = &union_zipped_chrom_interval_maps[&chrom];

                let binned_values: Vec<(I64Interval, Vec<Option<Value>>)> =
                    interval_maps
                        .iter()
                        .skip(1)
                        .fold(
                            interval_maps
                                .first()
                                .expect("interval maps cannot be empty")
                                .iter()
                                .into_binned_interval_iter(
                                    bin_size,
                                    AggregateOp::Average,
                                    Box::new(|item| (*item.0, *item.1)),
                                )
                                .into_common_refinement_zipped(),
                            |common_refinement, map| {
                                common_refinement.common_refinement_flat_zip(
                                    map.iter().into_binned_interval_iter(
                                        bin_size,
                                        AggregateOp::Average,
                                        Box::new(|item| (*item.0, *item.1)),
                                    ),
                                )
                            },
                        )
                        .collect();
                (chrom, binned_values)
            })
            .collect())
    }
}

impl LinearTrackMixer {
    pub fn write_to_bed_file(
        &self,
        path: &str,
    ) -> Result<(), biofile::error::Error> {
        let mut writer = BedWriter::new(path)?;

        for chrom in crate::util::get_sorted_keys(&self.chrom_to_binned_values)
            .into_iter()
        {
            let iter = LinearTrackMixerIter {
                binned_values_iter: self.chrom_to_binned_values[&chrom].iter(),
                weights: &self.weights,
            };
            let mut bed_data_line = BedDataLine {
                chrom,
                start: 0,
                end: 0,
                name: None,
                score: None,
                strand: None,
            };
            for (interval, value) in iter {
                bed_data_line.score = Some(value);
                bed_data_line.start = interval.get_start();
                // the end coordinate is exclusive under the BED format
                bed_data_line.end = interval.get_end() + 1;
                writer.write_bed_line(&bed_data_line)?;
            }
        }
        Ok(())
    }
}

impl<'a> Iterator for LinearTrackMixerIter<'a> {
    type Item = (I64Interval, Value);

    fn next(&mut self) -> Option<Self::Item> {
        match self.binned_values_iter.next() {
            None => None,
            Some((interval, values)) => {
                let value = self
                    .weights
                    .iter()
                    .zip(values.iter())
                    .fold(0., |acc, (w, v)| acc + w * v.unwrap_or(0.));
                Some((*interval, value))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::linear_track_mixer::{LinearTrackMixerIter, TrackZipper};
    use biofile::bed::Bed;
    use math::interval::I64Interval;
    use std::{
        collections::HashSet,
        io::{BufWriter, Write},
    };
    use tempfile::NamedTempFile;

    #[test]
    fn test_linear_track_mixer_iter() {
        let binned_values = vec![
            (I64Interval::new(0, 24), vec![Some(5.), None, Some(10.)]),
            (I64Interval::new(20, 49), vec![None, Some(20.), Some(100.)]),
            (I64Interval::new(50, 99), vec![None, Some(10.), None]),
            (I64Interval::new(0, 9), vec![Some(5.), Some(8.), Some(10.)]),
            (I64Interval::new(10, 19), vec![None, None, None]),
        ];
        let weights = vec![0.2, 0.5, 0.3];
        let mut iter = LinearTrackMixerIter {
            binned_values_iter: binned_values.iter(),
            weights: &weights,
        };
        assert_eq!(iter.next(), Some((I64Interval::new(0, 24), 4.)));
        assert_eq!(iter.next(), Some((I64Interval::new(20, 49), 40.)));
        assert_eq!(iter.next(), Some((I64Interval::new(50, 99), 5.)));
        assert_eq!(iter.next(), Some((I64Interval::new(0, 9), 8.)));
        assert_eq!(iter.next(), Some((I64Interval::new(10, 19), 0.)));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_linear_track_mixer() {
        let bed_1_path = {
            let bed_1 = NamedTempFile::new().unwrap();
            {
                let mut writer = BufWriter::new(&bed_1);
                writer
                    .write_fmt(format_args!(
                        "chr1 100 200 name_1 10\n\
                        chr1 150 250 name_2 5\n\
                        chr3 2500 2600 name_6 25\n"
                    ))
                    .unwrap();
            }
            bed_1.into_temp_path()
        };

        let bed_2_path = {
            let bed_2 = NamedTempFile::new().unwrap();
            {
                let mut writer = BufWriter::new(&bed_2);
                writer
                    .write_fmt(format_args!(
                        "chr1 100 200 name_1 100\n\
                        chr3 2500 2575 name_6 1000\n"
                    ))
                    .unwrap();
            }
            bed_2.into_temp_path()
        };

        let zipper = TrackZipper::new(
            vec![
                (0.4, Bed::new(bed_1_path.to_str().unwrap())),
                (0.6, Bed::new(bed_2_path.to_str().unwrap())),
            ],
            None,
        )
        .unwrap();

        let target_chroms: HashSet<String> =
            vec!["chr1".to_string(), "chr3".to_string()]
                .into_iter()
                .collect();

        let mixer = zipper.to_linear_track_mixer(&target_chroms, 50).unwrap();
        let mixed_bed_path = NamedTempFile::new().unwrap().into_temp_path();
        mixer
            .write_to_bed_file(mixed_bed_path.to_str().unwrap())
            .unwrap();
        let mixed_bed = Bed::new(mixed_bed_path.to_str().unwrap());
        let x = mixed_bed
            .get_chrom_to_interval_to_val::<f64, _>(None)
            .unwrap();

        {
            let chr1_map = &x["chr1"];
            let mut chr1_map_iter = chr1_map.iter();
            assert_eq!(
                chr1_map_iter.next(),
                Some((&I64Interval::new(100, 149), &64.))
            );
            assert_eq!(
                chr1_map_iter.next(),
                Some((&I64Interval::new(150, 199), &66.))
            );
            assert_eq!(
                chr1_map_iter.next(),
                Some((&I64Interval::new(200, 249), &2.))
            );
            assert_eq!(chr1_map_iter.next(), None);
        }
        {
            let chr3_map = &x["chr3"];
            let mut chr3_map_iter = chr3_map.iter();
            assert_eq!(
                chr3_map_iter.next(),
                Some((&I64Interval::new(2500, 2549), &610.))
            );
            assert_eq!(
                chr3_map_iter.next(),
                Some((&I64Interval::new(2550, 2599), &310.))
            );
            assert_eq!(chr3_map_iter.next(), None);
        }
    }
}
