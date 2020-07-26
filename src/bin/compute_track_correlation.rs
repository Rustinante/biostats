use biofile::bed::Bed;
use clap::{clap_app, Arg};
use math::{
    interval::I64Interval,
    iter::{
        binned_interval_iter::{AggregateOp, IntoBinnedIntervalIter},
        CommonRefinementZip, UnionZip,
    },
    set::traits::Finite,
    stats::correlation::weighted_correlation,
};
use program_flow::{
    argparse::{extract_optional_numeric_arg, extract_str_arg},
    OrExit,
};

fn main() {
    let mut app = clap_app!(compute_track_correlation =>
        (version: "0.1")
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

    let chrom_interval_map_a = bed_a
        .get_chrom_to_interval_to_val::<f64, _>()
        .unwrap_or_exit(Some(format_args!(
            "failed to get chrom interval map for {}",
            first_track_filepath
        )));
    let chrom_interval_map_b = bed_b
        .get_chrom_to_interval_to_val::<f64, _>()
        .unwrap_or_exit(Some(format_args!(
            "failed to get chrom interval map for {}",
            second_track_filepath
        )));

    for (chrom, map_list) in chrom_interval_map_a
        .union_zip(&chrom_interval_map_b)
        .into_iter()
    {
        let interval_map_a = map_list[0].unwrap();
        let interval_map_b = map_list[1].unwrap();
        let vec: Vec<(I64Interval, Vec<Option<f64>>)> = match bin_size {
            None => interval_map_a
                .iter()
                .common_refinement_zip(interval_map_b.iter())
                .collect(),
            Some(bin_size) => interval_map_a
                .iter()
                .into_binned_interval_iter(
                    bin_size,
                    AggregateOp::Sum,
                    Box::new(|item| (*item.0, *item.1)),
                )
                .common_refinement_zip(interval_map_b.iter().into_binned_interval_iter(
                    bin_size,
                    AggregateOp::Sum,
                    Box::new(|item| (*item.0, *item.1)),
                ))
                .collect(),
        };

        let correlation = weighted_correlation(
            || vec.iter(),
            |(interval, v)| {
                (
                    v[0].unwrap_or(0.),
                    v[1].unwrap_or(0.),
                    1. / interval.size() as f64,
                )
            },
        );
        println!("Pearson correlation: {:.5}", correlation);
    }
}
