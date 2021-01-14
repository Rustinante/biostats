use biofile::{
    bed::{Bed, BedDataLine, BedDataLineIter, BedWriter, Chrom},
    bedgraph::BedGraphDataLine,
};
use math::{
    interval::{traits::Interval, I64Interval},
    iter::{AggregateOp, IntoBinnedIntervalIter, WeightedSum},
    partition::integer_interval_map::IntegerIntervalMap,
    set::traits::{Intersect, Set},
    traits::ToIterator,
};
use num::{Float, FromPrimitive};
use std::{
    collections::{HashMap, HashSet},
    fmt::Debug,
    str::FromStr,
};

pub struct BedRefinery<D> {
    chrom_to_interval_map: HashMap<Chrom, IntegerIntervalMap<D>>,
    stats: RefineryStats,
}

pub struct RefineryStats {
    pub num_duplicate_lines: Option<i64>,
}

impl<D, E> BedRefinery<D>
where
    D: Float + FromPrimitive + FromStr<Err = E> + std::fmt::Display,
    E: Debug,
{
    pub fn new(
        track_filepath: &str,
        unique: bool,
        binarize_score: bool,
        filter_chroms: Option<HashSet<String>>,
        exclude_track_filepath: Option<String>,
        debug: bool,
    ) -> BedRefinery<D> {
        let exclude = if let Some(path) = exclude_track_filepath {
            // binarize_score is irrelevant for getting the intervals
            Some(Bed::new(&path, false).get_chrom_to_intervals())
        } else {
            None
        };

        let mut visited = HashSet::new();
        let mut num_pcr_duplicates = 0i64;

        let mut chrom_to_interval_map =
            HashMap::<Chrom, IntegerIntervalMap<D>>::new();

        let bed = Bed::new(track_filepath, binarize_score);
        for BedDataLine {
            chrom,
            start,
            end,
            name: _,
            score,
            strand,
        } in bed.to_iter(): BedDataLineIter<D>
        {
            if filter_chroms.is_some()
                && !filter_chroms.as_ref().unwrap().contains(&chrom)
            {
                continue;
            }

            let interval = I64Interval::new(start, end - 1);

            if let Some(chrom_to_excluded_intervals) = exclude.as_ref() {
                if let Some(excluded_intervals) =
                    chrom_to_excluded_intervals.get(&chrom)
                {
                    if interval
                        .has_non_empty_intersection_with(excluded_intervals)
                    {
                        continue;
                    }
                }
            }
            if unique {
                if !visited.insert((chrom.clone(), start, end, strand)) {
                    // duplicate PCR reads
                    num_pcr_duplicates += 1;
                    if debug {
                        eprintln!(
                            "PCR duplciate (chrom, start, end, strand): \
                            ({}, {}, {}, {:?})",
                            chrom, start, end, strand
                        )
                    }
                    continue;
                }
            }
            assert!(
                end > 0,
                "the end coordinate must be positive, encountered \
                (chrom, start, end): ({}, {}, {})",
                chrom,
                start,
                end,
            );

            let interval_map = chrom_to_interval_map
                .entry(chrom)
                .or_insert_with(IntegerIntervalMap::new);

            interval_map.aggregate(interval, score.unwrap_or(D::zero()));
        }
        BedRefinery {
            chrom_to_interval_map,
            stats: RefineryStats {
                num_duplicate_lines: if unique {
                    Some(num_pcr_duplicates)
                } else {
                    None
                },
            },
        }
    }

    pub fn write_refined_bed(
        &self,
        out_path: &str,
        bin_size: i64,
        normalize: bool,
        scaling: Option<D>,
        out_bedgraph: bool,
    ) -> Result<(), biofile::error::Error> {
        let mut writer = BedWriter::new(out_path)?;
        for chrom in crate::util::get_sorted_keys(&self.chrom_to_interval_map) {
            let interval_map = &self.chrom_to_interval_map[&chrom];

            macro_rules! get_interval_value_iter {
                () => {
                    if bin_size == 0 {
                        Box::new(
                            interval_map
                                .iter()
                                .map(|(&interval, &val)| (interval, val)),
                        ) as Box<dyn Iterator<Item = (I64Interval, D)>>
                    } else {
                        Box::new(interval_map.iter().into_binned_interval_iter(
                            bin_size,
                            AggregateOp::Average,
                            Box::new(|item| (*item.0, *item.1)),
                        )) as Box<dyn Iterator<Item = (I64Interval, D)>>
                    }
                };
            }

            let normalization_constant = if normalize {
                get_interval_value_iter!().weighted_sum()
            } else {
                D::one()
            };
            if normalization_constant == D::zero() {
                return Err(biofile::error::Error::Generic(
                    "cannot normalize the values when they sum to zero.".into(),
                ));
            }
            let scaling = scaling.unwrap_or(D::one()) / normalization_constant;

            if out_bedgraph {
                let mut bedgraph_line = BedGraphDataLine {
                    chrom: chrom.to_string(),
                    start: 0,
                    end_exclusive: 0,
                    value: D::zero(),
                };
                get_interval_value_iter!().try_for_each(
                    |(interval, value): (I64Interval, D)|
                        -> Result<(), biofile::error::Error>{
                        if !interval.is_empty() {
                            bedgraph_line.start = interval.get_start();
                            bedgraph_line.end_exclusive =
                                interval.get_end() + 1i64;

                            bedgraph_line.value = value * scaling;
                            writer.write_bedgraph_line(&bedgraph_line)?;
                        }
                        Ok(())
                    })?;
            } else {
                let mut bed_line = BedDataLine {
                    chrom: chrom.to_string(),
                    start: 0,
                    end: 0,
                    name: None,
                    score: None::<D>,
                    strand: None,
                };
                get_interval_value_iter!().try_for_each(
                    |(interval, value): (I64Interval, D)|
                        -> Result<(), biofile::error::Error>{
                        if !interval.is_empty() {
                            bed_line.start = interval.get_start();
                            // the end is exclusive in the BED format
                            bed_line.end = interval.get_end() + 1i64;
                            bed_line.score = Some(value * scaling);
                            writer.write_bed_line(&bed_line)?;
                        }
                        Ok(())
                    })?;
            }
        }
        Ok(())
    }

    pub fn get_chrom_to_interval_map(
        &self,
    ) -> &HashMap<Chrom, IntegerIntervalMap<D>> {
        &self.chrom_to_interval_map
    }

    pub fn stats(&self) -> &RefineryStats {
        &self.stats
    }
}

#[cfg(test)]
mod tests {
    use crate::{bed_refinery::BedRefinery, util::manifest_path_join};
    use math::interval::I64Interval;

    #[test]
    fn test_refinery() {
        macro_rules! interval_val {
            ($start:expr, $end:expr, $val:expr) => {
                (I64Interval::new($start, $end), $val)
            };
        }
        macro_rules! check_intervals {
            ($interval_map:expr, $expected_pairs:expr) => {
                assert_eq!($interval_map.len(), $expected_pairs.len());

                for ((interval, value), expected) in
                    $interval_map.iter().zip($expected_pairs.iter())
                {
                    assert_eq!(*interval, expected.0);
                    assert_almost_eq!(*value, expected.1);
                }
            };
        }
        {
            let refinery = BedRefinery::<f64>::new(
                manifest_path_join("tests/test_3.bed").to_str().unwrap(),
                false,
                false,
                None,
                false,
            );
            let chrom_to_interval_map = refinery.get_chrom_to_interval_map();
            let expected_chr1 = vec![
                interval_val!(0, 2, 3.),
                interval_val!(3, 4, 5.),
                interval_val!(5, 9, 3.),
                interval_val!(10, 19, 7.),
                interval_val!(20, 25, 20.),
                interval_val!(26, 26, 13.),
                interval_val!(34, 36, 9.),
            ];
            let chr1_interval_map = &chrom_to_interval_map["chr1"];
            check_intervals!(chr1_interval_map, expected_chr1);
        }

        {
            let refinery = BedRefinery::<f32>::new(
                manifest_path_join("tests/test_4.bed").to_str().unwrap(),
                false,
                false,
                None,
                false,
            );
            let chrom_to_interval_map = refinery.get_chrom_to_interval_map();
            let expected_chr1 = vec![
                interval_val!(2, 2, 1.),
                interval_val!(3, 3, 8.),
                interval_val!(4, 4, 7.),
                interval_val!(5, 7, 9.),
                interval_val!(8, 9, 2.),
                interval_val!(10, 20, 12.),
                interval_val!(21, 39, 2.),
            ];
            let chr1_interval_map = &chrom_to_interval_map["chr1"];
            check_intervals!(chr1_interval_map, expected_chr1);
        }
    }
}
