use biofile::bed::{Bed, BedDataLine, BedWriter, Chrom};
use math::{
    interval::{traits::Interval, I64Interval},
    iter::{AggregateOp, CommonRefinementZip, IntoBinnedIntervalIter},
};
use std::collections::{HashMap, HashSet};

type Coefficient = f64;
type Value = f64;

pub struct LinearTrackMixture {
    content: HashMap<Chrom, Vec<(I64Interval, Value)>>,
}

impl LinearTrackMixture {
    pub fn create(
        weighted_paths: Vec<(Coefficient, String)>,
        bin_size: i64,
        use_binary_score: bool,
        exclude_track_filepath: Option<String>,
        target_chroms: Option<HashSet<String>>,
    ) -> Result<LinearTrackMixture, biofile::error::Error> {
        let exclude = if let Some(path) = exclude_track_filepath {
            // binarize_score is irrelevant for getting the intervals
            Some(Bed::new(&path, false).get_chrom_to_intervals())
        } else {
            None
        };

        let weights: Vec<Coefficient> =
            weighted_paths.iter().map(|(w, _)| *w).collect();

        let (first_weight, rest_weights): (Coefficient, Vec<Coefficient>) = {
            let (first, rest) = weights
                .split_first()
                .expect("weighted_path cannot be empty");
            (*first, rest.iter().cloned().collect())
        };

        let beds: Vec<Bed> = weighted_paths
            .iter()
            .map(|(_, p)| Bed::new(p, use_binary_score))
            .collect();

        let init: HashMap<Chrom, Vec<(I64Interval, Value)>> = beds
            .first()
            .expect("weighted_paths cannot be empty")
            .get_chrom_to_interval_to_val::<Value, _>(exclude.as_ref())?
            .into_iter()
            .filter_map(|(chrom, interval_to_val)| {
                if target_chroms.is_none()
                    || target_chroms.as_ref().unwrap().contains(&chrom)
                {
                    let binned_intervals = interval_to_val
                        .iter()
                        .into_binned_interval_iter(
                            bin_size,
                            AggregateOp::Average,
                            Box::new(|item| (*item.0, *item.1)),
                        )
                        .into_iter()
                        .map(|(interval, value)| {
                            (interval, value * first_weight)
                        })
                        .collect::<Vec<(I64Interval, Value)>>();
                    Some((chrom, binned_intervals))
                } else {
                    None
                }
            })
            .collect();

        let content = beds.iter().skip(1).enumerate().try_fold(
            init,
            |acc_chrom_to_binned_interval_values, (i, bed)| {
                bed.get_chrom_to_interval_to_val::<Value, _>(exclude.as_ref())?
                    .into_iter()
                    .filter_map(|(chrom, interval_to_val)| {
                        if target_chroms.is_some()
                            && !target_chroms.as_ref().unwrap().contains(&chrom)
                        {
                            return None;
                        }
                        let interval_values: Vec<(I64Interval, Value)> =
                            interval_to_val
                                .iter()
                                .into_binned_interval_iter(
                                    bin_size,
                                    AggregateOp::Average,
                                    Box::new(|item| (*item.0, *item.1)),
                                )
                                .collect();

                        let w: Coefficient = rest_weights[i];
                        match acc_chrom_to_binned_interval_values.get(&chrom) {
                            None => Some(Ok((chrom, interval_values))),

                            Some(acc_binned_interval_values) => {
                                let interval_values =
                                    acc_binned_interval_values
                                        .iter()
                                        .into_binned_interval_iter(
                                            bin_size,
                                            AggregateOp::Average,
                                            Box::new(|&item| (item.0, item.1)),
                                        )
                                        .common_refinement_zip(
                                            interval_values
                                                .iter()
                                                .into_binned_interval_iter(
                                                    bin_size,
                                                    AggregateOp::Average,
                                                    Box::new(|&item| {
                                                        (item.0, item.1)
                                                    }),
                                                ),
                                        )
                                        .map(|(interval, values)| {
                                            let acc = values[0].unwrap_or(0f64);
                                            let val = values[1].unwrap_or(0f64);
                                            (interval, acc + val * w)
                                        })
                                        .collect();
                                Some(Ok((chrom, interval_values)))
                            }
                        }
                    })
                    .collect::<Result<
                        HashMap<Chrom, Vec<(I64Interval, Value)>>,
                        biofile::error::Error,
                    >>()
            },
        )?;
        Ok(LinearTrackMixture {
            content,
        })
    }

    pub fn write_to_bed_file(
        &self,
        path: &str,
    ) -> Result<(), biofile::error::Error> {
        let mut writer = BedWriter::new(path)?;

        for chrom in crate::util::get_sorted_keys(&self.content).into_iter() {
            let mut bed_data_line = BedDataLine {
                chrom: chrom.clone(),
                start: 0,
                end: 0,
                name: None,
                score: None,
                strand: None,
            };
            for (interval, value) in self.content[&chrom].iter() {
                bed_data_line.start = interval.get_start();
                // the end coordinate is exclusive under the BED format
                bed_data_line.end = interval.get_end() + 1;
                bed_data_line.score = Some(*value);
                writer.write_bed_line(&bed_data_line)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::linear_track_mixture::LinearTrackMixture;
    use biofile::bed::Bed;
    use math::interval::I64Interval;
    use std::{
        collections::HashSet,
        io::{BufWriter, Write},
    };
    use tempfile::NamedTempFile;

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

        macro_rules! check_chrom {
            (
                $iter:expr, $($start_end_val:expr), *
            ) => {
                {
                    $(
                        let start = $start_end_val.0;
                        let end = $start_end_val.1;
                        let val = $start_end_val.2;
                        assert_eq!(
                            $iter.next(),
                            Some((&I64Interval::new(start, end), &val))
                        );
                    )*
                }
                assert_eq!(
                    $iter.next(),
                    None
                );
            };
        }

        let bin_size = 50;

        {
            let mixture = LinearTrackMixture::create(
                vec![
                    (0.4, bed_1_path.to_str().unwrap().to_string()),
                    (0.6, bed_2_path.to_str().unwrap().to_string()),
                ],
                bin_size,
                false,
                None,
                Some(
                    vec!["chr1".to_string(), "chr3".to_string()]
                        .into_iter()
                        .collect::<HashSet<String>>(),
                ),
            )
            .unwrap();

            let mixed_path = NamedTempFile::new().unwrap().into_temp_path();
            mixture
                .write_to_bed_file(mixed_path.to_str().unwrap())
                .unwrap();
            let x = Bed::new(mixed_path.to_str().unwrap(), false)
                .get_chrom_to_interval_to_val::<f64, _>(None)
                .unwrap();

            let mut chr1_map_iter = x["chr1"].iter();
            check_chrom!(
                chr1_map_iter,
                (100, 149, 64.),
                (150, 199, 66.),
                (200, 249, 2.)
            );

            let mut chr3_map_iter = x["chr3"].iter();
            check_chrom!(chr3_map_iter, (2500, 2549, 610.), (2550, 2599, 310.));
        }

        {
            let mixture = LinearTrackMixture::create(
                vec![
                    (0.5, bed_1_path.to_str().unwrap().to_string()),
                    (0.5, bed_2_path.to_str().unwrap().to_string()),
                ],
                bin_size,
                false,
                None,
                None,
            )
            .unwrap();

            let mixed_path = NamedTempFile::new().unwrap().into_temp_path();
            mixture
                .write_to_bed_file(mixed_path.to_str().unwrap())
                .unwrap();
            let x = Bed::new(mixed_path.to_str().unwrap(), false)
                .get_chrom_to_interval_to_val::<f64, _>(None)
                .unwrap();

            let mut chr1_map_iter = x["chr1"].iter();
            check_chrom!(
                chr1_map_iter,
                (100, 149, 55.),
                (150, 199, 57.5),
                (200, 249, 2.5)
            );

            let mut chr3_map_iter = x["chr3"].iter();
            check_chrom!(
                chr3_map_iter,
                (2500, 2549, 512.5),
                (2550, 2599, 262.5)
            );
        }

        {
            let mixture = LinearTrackMixture::create(
                vec![
                    (10., bed_1_path.to_str().unwrap().to_string()),
                    (0.1, bed_2_path.to_str().unwrap().to_string()),
                ],
                bin_size,
                false,
                None,
                Some(
                    vec!["chr3".to_string()]
                        .into_iter()
                        .collect::<HashSet<String>>(),
                ),
            )
            .unwrap();

            let mixed_path = NamedTempFile::new().unwrap().into_temp_path();
            mixture
                .write_to_bed_file(mixed_path.to_str().unwrap())
                .unwrap();
            let x = Bed::new(mixed_path.to_str().unwrap(), false)
                .get_chrom_to_interval_to_val::<f64, _>(None)
                .unwrap();

            assert!(!x.contains_key("chr1"));
            assert!(x.contains_key("chr3"));

            let mut chr3_map_iter = x["chr3"].iter();
            check_chrom!(chr3_map_iter, (2500, 2549, 350.), (2550, 2599, 300.));
        }
    }
}
