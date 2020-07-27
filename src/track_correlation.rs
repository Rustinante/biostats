use biofile::bed::Bed;
use math::{
    interval::I64Interval,
    iter::{
        AggregateOp, BinnedIntervalIter, CommonRefinementZip, CommonRefinementZipped,
        ConcatenatedIter, IntoBinnedIntervalIter, UnionZip,
    },
    partition::integer_interval_map::IntegerIntervalMap,
    set::traits::Finite,
    stats::correlation::weighted_correlation,
};
use std::collections::HashSet;

/// The weight is the reciprocal of the interval size so as to produce the mean of the values in
/// the interval. Each bin is considered a single entity of "weight" 1 when taking correlations.
macro_rules! binned_extractor {
    () => {
        |(interval, v)| {
            (
                v[0].unwrap_or(0.),
                v[1].unwrap_or(0.),
                1. / interval.size() as f64,
            )
        }
    };
}

/// The weight is the size of the interval since each value has to be multiplied by the number of
/// elements in the interval as the proper weighting in the correlation calculation.
macro_rules! non_binned_extractor {
    () => {
        |(interval, v)| {
            (
                v[0].unwrap_or(0.),
                v[1].unwrap_or(0.),
                interval.size() as f64,
            )
        }
    };
}

pub type ChromCorrelations = Vec<(String, Vec<f64>)>;
pub type OverallCorrelations = Vec<f64>;

pub fn compute_track_correlations(
    first_track_filepath: &str,
    second_track_filepath: &str,
    bin_sizes: Vec<i64>,
    target_chroms: HashSet<String>,
) -> Result<(ChromCorrelations, OverallCorrelations), String> {
    let bed_a = Bed::new(&first_track_filepath);
    let bed_b = Bed::new(&second_track_filepath);

    eprintln!(
        "=> Constructing chrom interval map for {}",
        first_track_filepath
    );
    let chrom_interval_map_a = match bed_a.get_chrom_to_interval_to_val::<f64, _>() {
        Ok(bed) => bed,
        Err(why) => {
            return Err(format!(
                "failed to get chrom interval map for {}: {}",
                first_track_filepath, why
            ))
        }
    };

    eprintln!(
        "=> Constructing chrom interval map for {}",
        second_track_filepath
    );
    let chrom_interval_map_b = match bed_b.get_chrom_to_interval_to_val::<f64, _>() {
        Ok(bed) => bed,
        Err(why) => {
            return Err(format!(
                "failed to get chrom interval map for {}: {}",
                second_track_filepath, why
            ))
        }
    };

    let get_target_interval_maps = || {
        chrom_interval_map_a
            .union_zip(&chrom_interval_map_b)
            .into_iter()
            .filter_map(|(chrom, map_list)| {
                if target_chroms.contains(&chrom) {
                    let interval_map_a = map_list[0].unwrap();
                    let interval_map_b = map_list[1].unwrap();
                    Some((chrom, interval_map_a, interval_map_b))
                } else {
                    None
                }
            })
    };

    eprintln!("=> Computing correlations");
    let chrom_correlations: Vec<(String, Vec<f64>)> = get_target_interval_maps()
        .map(|(chrom, map_a, map_b)| {
            eprintln!("chrom: {}", chrom);
            let correlations = bin_sizes
                .iter()
                .map(|&s| match s {
                    0 => {
                        let vec: Vec<(I64Interval, Vec<Option<f64>>)> =
                            a_common_refine_b(map_a, map_b).collect();
                        weighted_correlation(|| vec.iter(), non_binned_extractor!())
                    }
                    non_zero => {
                        let vec: Vec<(I64Interval, Vec<Option<f64>>)> =
                            a_bin_b(map_a, map_b, non_zero).collect();
                        weighted_correlation(|| vec.iter(), binned_extractor!())
                    }
                })
                .collect();
            (chrom, correlations)
        })
        .collect();

    eprintln!("=> Computing overall correlations");
    let overall_correlations: Vec<f64> = bin_sizes
        .iter()
        .map(|&s| match s {
            0 => weighted_correlation(
                || {
                    ConcatenatedIter::from_iters(
                        get_target_interval_maps()
                            .map(|(_, map_a, map_b)| a_common_refine_b(map_a, map_b))
                            .collect(),
                    )
                },
                non_binned_extractor!(),
            ),
            non_zero => weighted_correlation(
                || {
                    ConcatenatedIter::from_iters(
                        get_target_interval_maps()
                            .map(|(_, map_a, map_b)| a_bin_b(map_a, map_b, non_zero))
                            .collect(),
                    )
                },
                binned_extractor!(),
            ),
        })
        .collect();
    Ok((chrom_correlations, overall_correlations))
}

type Boundary = i64;
type Value = f64;
type BtreeMapIter<'a> = std::collections::btree_map::Iter<'a, I64Interval, f64>;
type BinnedIter<'a> =
    BinnedIntervalIter<std::collections::btree_map::Iter<'a, I64Interval, f64>, f64>;

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
            AggregateOp::Sum,
            Box::new(|item| (*item.0, *item.1)),
        )
        .common_refinement_zip(map_b.iter().into_binned_interval_iter(
            bin_size,
            AggregateOp::Sum,
            Box::new(|item| (*item.0, *item.1)),
        ))
}