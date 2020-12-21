use biofile::bed::{Bed, Chrom};
use math::{
    interval::{traits::Interval, I64Interval},
    iter::{AggregateOp, CommonRefinementZip, IntoBinnedIntervalIter},
    partition::integer_interval_map::IntegerIntervalMap,
};
use std::{
    collections::{HashMap, HashSet},
    fs::OpenOptions,
    io::{BufWriter, Write},
};

type Value = f64;

pub struct TrackZipper {
    pub bed_files: Vec<Bed>,
    list_of_chrom_interval_maps: Vec<HashMap<Chrom, IntegerIntervalMap<Value>>>,
}

impl TrackZipper {
    pub fn new(
        bed_files: Vec<Bed>,
        exclude_track_filepath: Option<&str>,
    ) -> Result<Self, biofile::error::Error> {
        let exclude = if let Some(path) = exclude_track_filepath {
            // binarize_score is irrelevant for getting the intervals
            Some(Bed::new(path, false).get_chrom_to_intervals())
        } else {
            None
        };

        let list_of_chrom_interval_maps: Vec<
            HashMap<Chrom, IntegerIntervalMap<Value>>,
        > = bed_files
            .iter()
            .map(|bed| bed.get_chrom_to_interval_to_val(exclude.as_ref()))
            .collect::<Result<
                Vec<HashMap<Chrom, IntegerIntervalMap<Value>>>,
                biofile::error::Error,
            >>()?;

        Ok(TrackZipper {
            bed_files,
            list_of_chrom_interval_maps,
        })
    }

    pub fn num_tracks(&self) -> usize {
        self.list_of_chrom_interval_maps.len()
    }

    pub fn chrom_to_binned_zipped_values(
        &self,
        target_chroms: Option<&HashSet<Chrom>>,
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

    pub fn write_concatenated_tracks(
        &self,
        target_chroms: Option<&HashSet<Chrom>>,
        bin_size: i64,
        out_path: &str,
    ) -> Result<(), biofile::error::Error> {
        let chrom_to_binned_zipped_values: HashMap<
            Chrom,
            Vec<(I64Interval, Vec<Option<Value>>)>,
        > = self.chrom_to_binned_zipped_values(target_chroms, bin_size)?;

        let chroms: Vec<String> = {
            let mut keys: Vec<String> =
                chrom_to_binned_zipped_values.keys().cloned().collect();
            keys.sort();
            keys
        };

        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(out_path)?;
        let mut writer = BufWriter::new(file);

        for c in chroms.into_iter() {
            for (interval, values) in &chrom_to_binned_zipped_values[&c] {
                // note that the end coordinate is exclusive in the BED format
                write!(
                    &mut writer,
                    "{} {} {}",
                    interval.get_start(),
                    interval.get_end() + 1,
                    c,
                )?;
                for v in values.iter() {
                    write!(&mut writer, " {}", v.unwrap_or(0.))?;
                }
                write!(&mut writer, "\n")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::track_zipper::TrackZipper;
    use biofile::bed::Bed;
    use std::{
        fs::OpenOptions,
        io::{BufRead, BufReader, BufWriter, Write},
    };
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_concatenated_tracks() {
        let bed_1_path = {
            let bed_1 = NamedTempFile::new().unwrap();
            {
                let mut writer = BufWriter::new(&bed_1);
                writer
                    .write_fmt(format_args!(
                        "chr1 100 110 name_1 17\n\
                        chr1 200 212 name_2 50\n\
                        chr3 1000 1025 name_6 250\n"
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
                        "chr1 100 175 name_1 15\n\
                        chr3 1025 1026 name_6 20\n"
                    ))
                    .unwrap();
            }
            bed_2.into_temp_path()
        };
        let bed_3_path = {
            let bed_3 = NamedTempFile::new().unwrap();
            {
                let mut writer = BufWriter::new(&bed_3);
                writer
                    .write_fmt(format_args!(
                        "chr2 10 20 name_1 10\n\
                        chr1 90 105 1000 7\n"
                    ))
                    .unwrap();
            }
            bed_3.into_temp_path()
        };

        let zipper = TrackZipper::new(
            vec![
                Bed::new(bed_1_path.to_str().unwrap(), false),
                Bed::new(bed_2_path.to_str().unwrap(), false),
                Bed::new(bed_3_path.to_str().unwrap(), false),
            ],
            None,
        )
        .unwrap();

        let out_path = NamedTempFile::new().unwrap().into_temp_path();
        zipper
            .write_concatenated_tracks(None, 25, out_path.to_str().unwrap())
            .unwrap();

        let reader = BufReader::new(
            OpenOptions::new()
                .read(true)
                .open(out_path.to_str().unwrap())
                .unwrap(),
        );
        let lines: Vec<(i32, i32, String, Vec<f64>)> = reader
            .lines()
            .map(|line| {
                let line = line.unwrap();
                let mut iter = line.split_whitespace();
                let start: i32 = iter.next().unwrap().parse().unwrap();
                let end: i32 = iter.next().unwrap().parse().unwrap();
                let chrom = iter.next().unwrap().to_string();
                let values: Vec<f64> =
                    iter.map(|val| val.parse::<f64>().unwrap()).collect();
                (start, end, chrom, values)
            })
            .collect();

        let expected: Vec<(i32, i32, String, Vec<f64>)> = vec![
            (75, 100, "chr1".into(), vec![0., 0., 2.8]),
            (100, 125, "chr1".into(), vec![6.8, 15., 1.4]),
            (125, 150, "chr1".into(), vec![0., 15., 0.]),
            (150, 175, "chr1".into(), vec![0., 15., 0.]),
            (200, 225, "chr1".into(), vec![24., 0., 0.]),
            (0, 25, "chr2".into(), vec![0., 0., 4.]),
            (1000, 1025, "chr3".into(), vec![250., 0., 0.]),
            (1025, 1050, "chr3".into(), vec![0., 0.8, 0.]),
        ];

        assert_eq!(lines.len(), expected.len());
        for (actual, expected) in lines.into_iter().zip(expected.into_iter()) {
            let (start, end, chrom, values) = actual;
            let (expected_start, expected_end, expected_chrom, expected_values) =
                expected;
            assert_eq!(start, expected_start);
            assert_eq!(end, expected_end);
            assert_eq!(chrom, expected_chrom);
            assert_vec_almost_eq!(values, expected_values);
        }
    }
}
