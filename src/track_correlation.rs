use crate::{top_k::get_top_k, util::get_chrom_interval_map};
use biofile::{bed::Bed, util::TrackVariant};
use math::{
    interval::I64Interval,
    iter::{
        AggregateOp, BinnedIntervalIter, CommonRefinementZip,
        CommonRefinementZipped, ConcatenatedIter, IntoBinnedIntervalIter,
        UnionZip,
    },
    partition::integer_interval_map::IntegerIntervalMap,
    set::traits::Finite,
    stats::correlation::weighted_correlation,
};
use std::collections::HashSet;

/// The weight is the reciprocal of the interval size so as to produce the mean
/// of the values in the interval. Each bin is considered a single entity of
/// "weight" 1 when taking correlations.
macro_rules! binned_extractor {
    ($transform: expr, $transform_type: expr) => {
        |(interval, v)| {
            (
                $transform(v[0].unwrap_or(0.), $transform_type),
                $transform(v[1].unwrap_or(0.), $transform_type),
                interval.size() as f64,
            )
        }
    };
}

/// The weight is the size of the interval since each value has to be multiplied
/// by the number of elements in the interval as the proper weighting in the
/// correlation calculation.
macro_rules! non_binned_extractor {
    ($transform: expr, $transform_type: expr) => {
        |(interval, v)| {
            (
                $transform(v[0].unwrap_or(0.), $transform_type),
                $transform(v[1].unwrap_or(0.), $transform_type),
                interval.size() as f64,
            )
        }
    };
}

pub type ChromCorrelations = Vec<(String, Vec<f64>)>;
pub type OverallCorrelations = Vec<f64>;

type Coord = i64;

pub fn compute_track_correlations(
    first_track: &TrackVariant,
    second_track: &TrackVariant,
    bin_sizes: &Vec<Coord>,
    target_chroms: Option<HashSet<String>>,
    value_transform: ValueTransform,
    top_k: Option<i64>,
    exclude_track_filepath: Option<String>,
) -> Result<(ChromCorrelations, OverallCorrelations), String> {
    let exclude = if let Some(path) = exclude_track_filepath {
        // binarize_score is irrelevant for getting the intervals
        Some(Bed::new(&path, false).get_chrom_to_intervals())
    } else {
        None
    };

    eprintln!("=> Constructing chrom interval map for the first track");
    let chrom_interval_map_a =
        get_chrom_interval_map(first_track, exclude.as_ref())?;

    eprintln!("=> Constructing chrom interval map for the second track");
    let chrom_interval_map_b =
        get_chrom_interval_map(second_track, exclude.as_ref())?;

    let empty_interval_map = IntegerIntervalMap::new();
    let get_target_interval_maps = || {
        chrom_interval_map_a
            .union_zip(&chrom_interval_map_b)
            .into_iter()
            .filter_map(|(chrom, map_list)| {
                if target_chroms.is_none()
                    || target_chroms.as_ref().unwrap().contains(&chrom)
                {
                    let interval_map_a =
                        map_list[0].unwrap_or_else(|| &empty_interval_map);
                    let interval_map_b =
                        map_list[1].unwrap_or_else(|| &empty_interval_map);
                    Some((chrom, interval_map_a, interval_map_b))
                } else {
                    None
                }
            })
    };

    let get_a_bin_b_zipped =
        |map_a,
         map_b,
         bin_size|
         -> Result<Vec<(I64Interval, Vec<Option<f64>>)>, String> {
            if let Some(k) = top_k {
                let map_a_top_k = get_top_k(map_a, k, bin_size)?;
                let map_b_top_k = get_top_k(map_b, k, bin_size)?;
                Ok(a_bin_b(&map_a_top_k, &map_b_top_k, bin_size).collect())
            } else {
                Ok(a_bin_b(map_a, map_b, bin_size).collect())
            }
        };

    let chrom_correlations: Vec<(String, Vec<f64>)> =
        get_target_interval_maps()
            .map(|(chrom, map_a, map_b)| {
                eprintln!("=> Computing correlations for {}", chrom);

                let correlations: Result<Vec<f64>, String> = bin_sizes
                    .iter()
                    .map(|&s| match s {
                        0 => {
                            let vec: Vec<(I64Interval, Vec<Option<f64>>)> =
                                a_common_refine_b(map_a, map_b).collect();

                            Ok(weighted_correlation(
                                || vec.iter(),
                                non_binned_extractor!(
                                    apply_transform,
                                    value_transform
                                ),
                            ))
                        }
                        non_zero => {
                            let vec: Vec<(I64Interval, Vec<Option<f64>>)> =
                                get_a_bin_b_zipped(map_a, map_b, non_zero)?;

                            Ok(weighted_correlation(
                                || vec.iter(),
                                binned_extractor!(
                                    apply_transform,
                                    value_transform
                                ),
                            ))
                        }
                    })
                    .collect();

                correlations.map(|c| (chrom, c))
            })
            .collect::<Result<Vec<(String, Vec<f64>)>, String>>()?;

    eprintln!("=> Computing overall correlations");
    let overall_correlations: Vec<f64> = bin_sizes
        .iter()
        .map(|&s| match s {
            0 => weighted_correlation(
                || {
                    ConcatenatedIter::from_iters(
                        get_target_interval_maps()
                            .map(|(_, map_a, map_b)| {
                                a_common_refine_b(map_a, map_b)
                            })
                            .collect(),
                    )
                },
                non_binned_extractor!(apply_transform, value_transform),
            ),
            non_zero => weighted_correlation(
                || {
                    ConcatenatedIter::from_iters(
                        get_target_interval_maps()
                            .map(|(_, map_a, map_b)| {
                                a_bin_b(map_a, map_b, non_zero)
                            })
                            .collect(),
                    )
                },
                binned_extractor!(apply_transform, value_transform),
            ),
        })
        .collect();
    Ok((chrom_correlations, overall_correlations))
}

/// `Idenitty` does not change the value.
/// `Log` transforms any value x into sign(x) * ln(|x| + 1)
/// `Thresholding(t)` will restrict the absolute value to less than or equal to
/// `t`. `LogThresholding(t)` will first apply the thresholding and then apply
/// the log transform.
#[derive(Copy, Clone, PartialEq)]
pub enum ValueTransform {
    Identity,
    Thresholding(f64),
}

fn apply_transform(value: f64, transform: ValueTransform) -> f64 {
    match transform {
        ValueTransform::Identity => value,
        ValueTransform::Thresholding(t) => {
            if value > t {
                t
            } else if value < -t {
                -t
            } else {
                value
            }
        }
    }
}

type Boundary = i64;
type Value = f64;
type BtreeMapIter<'a> = std::collections::btree_map::Iter<'a, I64Interval, f64>;
type BinnedIter<'a> = BinnedIntervalIter<
    std::collections::btree_map::Iter<'a, I64Interval, f64>,
    f64,
>;

fn a_common_refine_b<'a>(
    map_a: &'a IntegerIntervalMap<f64>,
    map_b: &'a IntegerIntervalMap<f64>,
) -> CommonRefinementZipped<
    Boundary,
    BtreeMapIter<'a>,
    <BtreeMapIter<'a> as Iterator>::Item,
    I64Interval,
    Value,
> {
    map_a.iter().common_refinement_zip(map_b.iter())
}

fn a_bin_b<'a>(
    map_a: &'a IntegerIntervalMap<f64>,
    map_b: &'a IntegerIntervalMap<f64>,
    bin_size: i64,
) -> CommonRefinementZipped<
    Boundary,
    BinnedIter<'a>,
    <BinnedIter<'a> as Iterator>::Item,
    I64Interval,
    Value,
> {
    map_a
        .iter()
        .into_binned_interval_iter(
            bin_size,
            AggregateOp::Average,
            Box::new(|item| (*item.0, *item.1)),
        )
        .common_refinement_zip(map_b.iter().into_binned_interval_iter(
            bin_size,
            AggregateOp::Average,
            Box::new(|item| (*item.0, *item.1)),
        ))
}
