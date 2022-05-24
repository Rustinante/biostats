use biofile::{bed::Bed, bedgraph::BedGraph, util::TrackVariant};
use biostats::{
    top_k_overlap::{
        get_top_k_fraction_overlap_ratio,
        get_top_k_fraction_overlap_ratio_across_chroms,
    },
    util::{get_chrom_interval_map, get_excluded_interval_maps},
};
use clap::{clap_app, Arg};
use math::{
    iter::UnionZip, partition::integer_interval_map::IntegerIntervalMap,
};
use program_flow::{
    argparse::{
        extract_boolean_flag, extract_numeric_arg, extract_optional_str_arg,
        extract_optional_str_vec_arg, extract_str_arg,
    },
    debug_eprint_named_vars, eprint_named_vars, OrExit,
};

const ZERO_BIN_SIZE_STR: &str = "0";
const BINARIZE_SCORE: bool = false;

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
            Arg::with_name("top_k_fraction")
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
    let top_k_fraction = extract_numeric_arg(&matches, "top_k_fraction")
        .unwrap_or_exit(Some("failed to parse top-k"));

    let default_human_chroms =
        extract_boolean_flag(&matches, "default_human_chroms");

    let exclude = extract_optional_str_arg(&matches, "exclude");

    let first_bedgraph = extract_boolean_flag(&matches, "first_bedgraph");
    let second_bedgraph = extract_boolean_flag(&matches, "second_bedgraph");

    let first_track: TrackVariant = if first_bedgraph {
        TrackVariant::BedGraph(BedGraph::new(
            &first_track_filepath,
            BINARIZE_SCORE,
        ))
    } else {
        TrackVariant::Bed(Bed::new(&first_track_filepath, BINARIZE_SCORE))
    };

    let second_track: TrackVariant = if second_bedgraph {
        TrackVariant::BedGraph(BedGraph::new(
            &second_track_filepath,
            BINARIZE_SCORE,
        ))
    } else {
        TrackVariant::Bed(Bed::new(&second_track_filepath, BINARIZE_SCORE))
    };

    eprint_named_vars!(
        first_track_filepath,
        second_track_filepath,
        default_human_chroms,
        first_bedgraph,
        second_bedgraph
    );
    debug_eprint_named_vars!(exclude, bin_sizes, chroms, top_k_fraction);

    let exlcude = get_excluded_interval_maps(exclude);

    let chrom_interval_map_1 =
        get_chrom_interval_map(&first_track, exlcude.as_ref())
            .unwrap_or_exit(None::<String>);

    let chrom_interval_map_2 =
        get_chrom_interval_map(&second_track, exlcude.as_ref())
            .unwrap_or_exit(None::<String>);

    let invalid_bin_sizes: Vec<i64> = bin_sizes
        .iter()
        .filter(|&&s| s <= 0i64)
        .map(|&s| s)
        .collect();

    if !invalid_bin_sizes.is_empty() {
        eprintln!(
            "bin sizes must be positive, but these are not: {:?}",
            invalid_bin_sizes
        );
        std::process::exit(1);
    }

    let empty_interval_map = IntegerIntervalMap::new();
    for b in bin_sizes {
        println!(
            "=> computing top {} overlap with bin size {}",
            top_k_fraction, b
        );
        chrom_interval_map_1
            .union_zip(&chrom_interval_map_2)
            .into_iter()
            .for_each(|(chrom, map_list)| {
                let ratio = get_top_k_fraction_overlap_ratio(
                    &map_list[0].unwrap_or_else(|| &empty_interval_map),
                    &map_list[1].unwrap_or_else(|| &empty_interval_map),
                    top_k_fraction,
                    b,
                )
                .unwrap_or_exit(None::<String>);

                println!("{}, {} ", chrom, ratio);
            });

        let ratio = get_top_k_fraction_overlap_ratio_across_chroms(
            &chrom_interval_map_1,
            &chrom_interval_map_2,
            top_k_fraction,
            b,
        )
        .unwrap_or_exit(Some("failed to compute overall overlap ratio"));

        println!(
            "=> computing overall {} overlap with bin size {}",
            top_k_fraction, b
        );
        println!("top_k_overall, {}", ratio);
    }
}
