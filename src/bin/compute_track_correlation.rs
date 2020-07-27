use biofile::bed::Bed;
use biostats::util::get_default_human_chrom_inclusion_set;
use clap::{clap_app, Arg};
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
use program_flow::{
    argparse::{extract_optional_numeric_arg, extract_str_arg},
    OrExit,
};

fn main() {
    let mut app = clap_app!(compute_track_correlation =>
        (about: "Computes Pearson correlation between two tracks in BED format")
    );
    app = app
        .arg(
            Arg::with_name("first_track_filepath")
                .long("first")
                .short("a")
                .takes_value(true)
                .help("filepath to the first track"),
        )
        .arg(
            Arg::with_name("second_track_filepath")
                .long("second")
                .short("b")
                .takes_value(true)
                .help("filepath to the second track"),
        )
        .arg(
            Arg::with_name("bin_size")
                .long("bin")
                .takes_value(true)
                .help("group the basepairs into consecutive bins of size bin_size"),
        );
    let matches = app.get_matches();
    let first_track_filepath = extract_str_arg(&matches, "first_track_filepath");
    let second_track_filepath = extract_str_arg(&matches, "second_track_filepath");
    let bin_size = extract_optional_numeric_arg::<i64>(&matches, "bin_size")
        .unwrap_or_exit(Some("failed to parse bin_size"));
    println!(
        "first_track_filepath: {}\n\
        second_track_filepath: {}",
        first_track_filepath, second_track_filepath,
    );
    let bed_a = Bed::new(&first_track_filepath).unwrap_or_exit(Some(format!(
        "failed to read the BED file {}",
        first_track_filepath
    )));
    let bed_b = Bed::new(&second_track_filepath).unwrap_or_exit(Some(format!(
        "failed to read the BED file {}",
        first_track_filepath
    )));

    println!(
        "=> Constructing chrom interval map for {}",
        first_track_filepath
    );
    let chrom_interval_map_a = bed_a
        .get_chrom_to_interval_to_val::<f64, _>()
        .unwrap_or_exit(Some(format_args!(
            "failed to get chrom interval map for {}",
            first_track_filepath
        )));

    println!(
        "=> Constructing chrom interval map for {}",
        second_track_filepath
    );
    let chrom_interval_map_b = bed_b
        .get_chrom_to_interval_to_val::<f64, _>()
        .unwrap_or_exit(Some(format_args!(
            "failed to get chrom interval map for {}",
            second_track_filepath
        )));

    let extractor = |(interval, v): &(I64Interval, Vec<Option<f64>>)| -> (f64, f64, f64) {
        (
            v[0].unwrap_or(0.),
            v[1].unwrap_or(0.),
            1. / interval.size() as f64,
        )
    };

    println!("=> Computing correlations");
    let target_chroms = get_default_human_chrom_inclusion_set();
    let chrom_correlations: Vec<(String, f64)> = chrom_interval_map_a
        .union_zip(&chrom_interval_map_b)
        .into_iter()
        .filter_map(|(chrom, map_list)| {
            if !target_chroms.contains(&chrom) {
                println!("- skipping chromosome {}", chrom);
                None
            } else {
                let interval_map_a = map_list[0].unwrap();
                let interval_map_b = map_list[1].unwrap();
                let vec: Vec<(I64Interval, Vec<Option<f64>>)> = match bin_size {
                    None => map_a_common_refine_map_b(interval_map_a, interval_map_b).collect(),
                    Some(bin_size) => {
                        map_a_bin_map_b(interval_map_a, interval_map_b, bin_size).collect()
                    }
                };
                Some((chrom, weighted_correlation(|| vec.iter(), extractor)))
            }
        })
        .collect();

    let overall_correlation = match bin_size {
        None => weighted_correlation(
            || {
                ConcatenatedIter::from_iters(
                    chrom_interval_map_a
                        .union_zip(&chrom_interval_map_b)
                        .into_iter()
                        .filter_map(|(chrom, map_list)| {
                            if target_chroms.contains(&chrom) {
                                let interval_map_a = map_list[0].unwrap();
                                let interval_map_b = map_list[1].unwrap();
                                Some(map_a_common_refine_map_b(interval_map_a, interval_map_b))
                            } else {
                                None
                            }
                        })
                        .collect(),
                )
            },
            |(interval, v)| {
                (
                    v[0].unwrap_or(0.),
                    v[1].unwrap_or(0.),
                    1. / interval.size() as f64,
                )
            },
        ),
        Some(bin_size) => weighted_correlation(
            || {
                ConcatenatedIter::from_iters(
                    chrom_interval_map_a
                        .union_zip(&chrom_interval_map_b)
                        .into_iter()
                        .filter_map(|(chrom, map_list)| {
                            if target_chroms.contains(&chrom) {
                                let interval_map_a = map_list[0].unwrap();
                                let interval_map_b = map_list[1].unwrap();
                                Some(map_a_bin_map_b(interval_map_a, interval_map_b, bin_size))
                            } else {
                                None
                            }
                        })
                        .collect(),
                )
            },
            |(interval, v)| {
                (
                    v[0].unwrap_or(0.),
                    v[1].unwrap_or(0.),
                    1. / interval.size() as f64,
                )
            },
        ),
    };

    for (chrom, correlation) in chrom_correlations.iter() {
        println!("{} Pearson correlation: {:.5}", chrom, correlation);
    }
    println!("overall Pearson correlation: {:.5}", overall_correlation);
}

type CommonBtreeRefinement<'a> = CommonRefinementZipped<
    i64,
    std::collections::btree_map::Iter<'a, I64Interval, f64>,
    (&'a I64Interval, &'a f64),
    I64Interval,
    f64,
>;

type CommonBinnedRefinement<'a> = CommonRefinementZipped<
    i64,
    BinnedIntervalIter<std::collections::btree_map::Iter<'a, I64Interval, f64>, f64>,
    (I64Interval, f64),
    I64Interval,
    f64,
>;

fn map_a_common_refine_map_b<'a>(
    map_a: &'a IntegerIntervalMap<f64>,
    map_b: &'a IntegerIntervalMap<f64>,
) -> CommonBtreeRefinement<'a> {
    map_a.iter().common_refinement_zip(map_b.iter())
}

fn map_a_bin_map_b<'a>(
    map_a: &'a IntegerIntervalMap<f64>,
    map_b: &'a IntegerIntervalMap<f64>,
    bin_size: i64,
) -> CommonBinnedRefinement<'a> {
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
