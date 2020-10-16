use biostats::{
    track_correlation::{compute_track_correlations, ValueTransform},
    util::get_default_human_chrom_inclusion_set,
};
use clap::{clap_app, Arg};
use program_flow::{
    argparse::{
        extract_boolean_flag, extract_optional_numeric_arg,
        extract_optional_str_arg, extract_optional_str_vec_arg,
        extract_str_arg,
    },
    debug_eprint_named_vars, eprint_named_vars, OrExit,
};
use std::{collections::HashSet, iter::FromIterator};

const ZERO_BIN_SIZE_STR: &str = "0";

fn main() {
    let mut app = clap_app!(compute_track_correlation =>
        (about: "Computes the Pearson correlation between two tracks stored in \
        BED format")
    );
    app = app
        .arg(
            Arg::with_name("first_track_filepath")
                .long("first")
                .short("a")
                .takes_value(true)
                .help("filepath to the first track."),
        )
        .arg(
            Arg::with_name("second_track_filepath")
                .long("second")
                .short("b")
                .takes_value(true)
                .help("filepath to the second track."),
        )
        .arg(
            Arg::with_name("bin_size")
                .long("bin")
                .takes_value(true)
                .multiple(true)
                .long_help(
                    "Group the base pairs into consecutive bins of size \
                    bin_size, aligned at index 0. A bin size of 0 means not to \
                    bin. A correlation will be computed for each provided bin \
                    value.",
                ),
        )
        .arg(
            Arg::with_name("chroms")
                .long("chroms")
                .short("c")
                .takes_value(true)
                .multiple(true)
                .long_help(
                    "Only process the specified chromosome. Specify this \
                    multiple times if you want to restrict to multiple \
                    chromosomes, e.g., --chroms chr1 --chroms chr2",
                ),
        )
        .arg(
            Arg::with_name("default_human_chroms")
                .long("default-human-chrom")
                .short("d")
                .conflicts_with("chroms")
                .help(
                    "Only computes histograms for chromosome chr1, chr2, ... \
                    chr22, chrX, chrY.",
                ),
        )
        .arg(
            Arg::with_name("threshold")
                .long("threshold")
                .short("t")
                .takes_value(true)
                .help(
                    "Apply thresholding to the aggregate value x at each base \
                    pair, while preserving the sign, i.e., \
                    x => sign(x) * min(|x|, threshold)) \
                    Threshold must be positive.",
                ),
        )
        .arg(
            Arg::with_name("exclude")
                .long("exclude")
                .short("v")
                .takes_value(true)
                .help(
                    "Path to a BED-like file where only the chromosome, start \
                    and end fields are required. Lines from other BED files \
                    that overlap with any of the coordinates in this 'exclude' \
                    file will be ignroed when computing correlations.",
                ),
        );
    let matches = app.get_matches();
    let first_track_filepath =
        extract_str_arg(&matches, "first_track_filepath");
    let second_track_filepath =
        extract_str_arg(&matches, "second_track_filepath");
    let bin_sizes: Vec<i64> =
        extract_optional_str_vec_arg(&matches, "bin_size")
            .unwrap_or_else(|| vec![ZERO_BIN_SIZE_STR.to_string()])
            .into_iter()
            .map(|s| {
                s.parse::<i64>().unwrap_or_exit(Some(format_args!(
                    "failed to parse {}",
                    &s
                )))
            })
            .collect();
    let chroms = extract_optional_str_vec_arg(&matches, "chroms");

    let default_human_chroms =
        extract_boolean_flag(&matches, "default_human_chroms");

    let log_transform = extract_boolean_flag(&matches, "log_transform");
    let threshold = extract_optional_numeric_arg(&matches, "threshold")
        .unwrap_or_exit(Some("failed to parse threshold"));
    let exclude = extract_optional_str_arg(&matches, "exclude");
    eprint_named_vars!(
        first_track_filepath,
        second_track_filepath,
        default_human_chroms,
        log_transform
    );
    debug_eprint_named_vars!(threshold, exclude, bin_sizes, chroms);

    let target_chroms = match chroms {
        Some(chroms) => {
            if default_human_chroms {
                eprintln!(
                    "at most one of chroms and default_human_chroms can be specified"
                );
                std::process::exit(1);
            } else {
                Some(HashSet::from_iter(chroms.into_iter()))
            }
        }
        None => {
            if default_human_chroms {
                Some(get_default_human_chrom_inclusion_set())
            } else {
                None
            }
        }
    };

    let transform_type = if let Some(t) = threshold {
        ValueTransform::Thresholding(t)
    } else {
        ValueTransform::Identity
    };

    let (chrom_correlations, overall_correlations) =
        compute_track_correlations(
            &first_track_filepath,
            &second_track_filepath,
            bin_sizes,
            target_chroms,
            transform_type,
            exclude,
        )
        .unwrap_or_exit(Some("failed to compute track correlations"));

    for (chrom, correlation) in chrom_correlations.iter() {
        print!("{}, ", chrom);
        correlation.iter().for_each(|c| print!("{:.5}, ", c));
        println!();
    }
    print!("overall, ");
    overall_correlations
        .iter()
        .for_each(|c| print!("{:.5}, ", c));
    println!();
}
