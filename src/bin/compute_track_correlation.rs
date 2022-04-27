use biofile::{bed::Bed, bedgraph::BedGraph, util::TrackVariant};
use biostats::{
    top_k_overlap::get_top_k_overlap_ratio,
    track_correlation::{compute_track_correlations, ValueTransform},
    util::{
        get_chrom_interval_map, get_default_human_chrom_inclusion_set,
        get_excluded_interval_maps,
    },
};
use clap::{clap_app, Arg};
use math::{
    iter::UnionZip, partition::integer_interval_map::IntegerIntervalMap,
};
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
            Arg::with_name("top_k")
                .long("top-k")
                .takes_value(true)
                .long_help(
                    "For each chromosome, compute the overlap ratio \
                    between the top K bins in each track when \
                    the bin size is non-zero. When the bin size is 0 this \
                    argument has no effect. Currently does not apply to \
                    the overall correlations across chromosomes.",
                ),
        )
        .arg(Arg::with_name("binarize_score").long("binarize").help(
            "Each line in the original BED files will contribute a \
            unit score for the corresponding interval",
        ))
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
        )
        .arg(
            Arg::with_name("first_bedgraph")
                .long("first-bedgraph")
                .help(
                    "A flag to indicate that the first track is in the \
                    bedgraph format",
                ),
        )
        .arg(
            Arg::with_name("second_bedgraph")
                .long("second-bedgraph")
                .help(
                    "A flag to indicate that the second track is in the \
                    bedgraph format",
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
    let top_k = extract_optional_numeric_arg(&matches, "top_k")
        .unwrap_or_exit(Some("failed to parse top-k"));

    let binarize_score = extract_boolean_flag(&matches, "binarize_score");

    let default_human_chroms =
        extract_boolean_flag(&matches, "default_human_chroms");

    let log_transform = extract_boolean_flag(&matches, "log_transform");
    let threshold = extract_optional_numeric_arg(&matches, "threshold")
        .unwrap_or_exit(Some("failed to parse threshold"));
    let exclude = extract_optional_str_arg(&matches, "exclude");

    let first_bedgraph = extract_boolean_flag(&matches, "first_bedgraph");
    let second_bedgraph = extract_boolean_flag(&matches, "second_bedgraph");

    let first_track: TrackVariant = if first_bedgraph {
        TrackVariant::BedGraph(BedGraph::new(
            &first_track_filepath,
            binarize_score,
        ))
    } else {
        TrackVariant::Bed(Bed::new(&first_track_filepath, binarize_score))
    };

    let second_track: TrackVariant = if second_bedgraph {
        TrackVariant::BedGraph(BedGraph::new(
            &second_track_filepath,
            binarize_score,
        ))
    } else {
        TrackVariant::Bed(Bed::new(&second_track_filepath, binarize_score))
    };

    eprint_named_vars!(
        first_track_filepath,
        second_track_filepath,
        binarize_score,
        default_human_chroms,
        log_transform,
        first_bedgraph,
        second_bedgraph
    );
    debug_eprint_named_vars!(threshold, exclude, bin_sizes, chroms, top_k);

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
            &first_track,
            &second_track,
            &bin_sizes,
            target_chroms,
            transform_type,
            None,
            exclude.clone(),
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

    if let Some(top_k) = top_k {
        let exlcude = get_excluded_interval_maps(exclude);

        let chrom_interval_map_1 =
            get_chrom_interval_map(&first_track, exlcude.as_ref())
                .unwrap_or_exit(None::<String>);

        let chrom_interval_map_2 =
            get_chrom_interval_map(&second_track, exlcude.as_ref())
                .unwrap_or_exit(None::<String>);

        let empty_interval_map = IntegerIntervalMap::new();

        for b in bin_sizes {
            if b <= 0 {
                continue;
            }
            println!("top {} overlap with bin size {}", top_k, b);
            chrom_interval_map_1
                .union_zip(&chrom_interval_map_2)
                .into_iter()
                .for_each(|(chrom, map_list)| {
                    let ratio = get_top_k_overlap_ratio(
                        &map_list[0].unwrap_or_else(|| &empty_interval_map),
                        &map_list[1].unwrap_or_else(|| &empty_interval_map),
                        top_k,
                        b,
                    )
                    .unwrap_or_exit(None::<String>);

                    println!("chrom {}: {} ", chrom, ratio);
                });
        }
    }
}
