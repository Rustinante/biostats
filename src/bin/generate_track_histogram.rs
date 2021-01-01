use biostats::{
    track_histogram::generate_track_histograms,
    util::{extract_chrom_names, get_default_human_chrom_inclusion_set},
};
use clap::{clap_app, Arg};
use program_flow::{
    argparse::{
        extract_boolean_flag, extract_numeric_arg,
        extract_optional_numeric_arg, extract_optional_str_arg,
        extract_str_arg,
    },
    debug_eprint_named_vars, eprint_named_vars, OrExit,
};

const DEFAULT_HISTOGRAM_NUM_BUCKETS: usize = 10;

fn main() {
    let mut app = clap_app!(generate_track_histogram =>
        (about: "Generates histograms for the track stored in BED format")
    );
    app = app
        .arg(
            Arg::with_name("track_filepath")
                .takes_value(true)
                .help("BED filepath for the track"),
        )
        .arg(
            Arg::with_name("bin_size")
                .long("bin")
                .short("b")
                .takes_value(true)
                .required(true)
                .long_help(
                    "Group the base pairs into consecutive bins of size bin_size, \
                    aligned at index 0. A bin size of 0 means not to bin.",
                ),
        )
        .arg(
            Arg::with_name("binarize_score")
                .long("binarize")
                .help(
                    "Each line in the original BED files will contribute a \
                    unit score for the corresponding interval",
                ),
        )
        .arg(
            Arg::with_name("min")
                .long("min")
                .takes_value(true)
                .required(true)
                .help("Minimum value of the histogram, defaults to 0."),
        )
        .arg(
            Arg::with_name("max")
                .long("max")
                .takes_value(true)
                .required(true)
                .help("Maximum value of the histogram, required."),
        )
        .arg(
            Arg::with_name("num_histogram_buckets")
                .long("buckets")
                .takes_value(true)
                .help("Number of histogram buckets, defaults to 10."),
        )
        .arg(
            Arg::with_name("default_human_chrom")
                .long("default-human-chrom")
                .short("d")
                .conflicts_with("filter_chrom")
                .help("Only computes histograms for chromosome chr1, chr2, ... chr22, chrX, chrY"),
        )
        .arg(
            Arg::with_name("filter_chrom")
                .long("filter-chrom")
                .short("f")
                .takes_value(true)
                .help(
                    "Path to a file in which each line is the name of a chromosome. If provided, \
                    will only compute histograms for the specified list of chromosome names. \
                    The leading and trailing whitespaces of each line will be ignored",
                ),
        );
    let matches = app.get_matches();
    let track_filepath = extract_str_arg(&matches, "track_filepath");
    let bin_size: i64 = extract_numeric_arg(&matches, "bin_size")
        .unwrap_or_exit(Some(format_args!("failed to parse --bin")));

    let binarize_score = extract_boolean_flag(&matches, "binarize_score");

    let min = extract_numeric_arg::<f64>(&matches, "min")
        .unwrap_or_exit(Some("failed to parse --min"));

    let max: f64 = extract_numeric_arg(&matches, "max")
        .unwrap_or_exit(Some("failed to parse --max"));

    let num_histogram_buckets = extract_optional_numeric_arg::<usize>(
        &matches,
        "num_histogram_buckets",
    )
    .unwrap_or_exit(Some("failed to parse --buckets"));
    let default_human_chrom =
        extract_boolean_flag(&matches, "default_human_chrom");
    let filter_chrom = extract_optional_str_arg(&matches, "filter_chrom");

    eprint_named_vars!(
        track_filepath,
        bin_size,
        binarize_score,
        min,
        max,
        default_human_chrom
    );
    debug_eprint_named_vars!(num_histogram_buckets, filter_chrom);

    let filter_chroms = if default_human_chrom {
        Some(get_default_human_chrom_inclusion_set())
    } else {
        if let Some(filter_chrom_filepath) = filter_chrom {
            println!("filter chromosomes:");
            Some(
                extract_chrom_names(&filter_chrom_filepath)
                    .unwrap_or_exit(Some(format_args!(
                        "failed to extract chromosome names from file {}",
                        filter_chrom_filepath
                    )))
                    .into_iter()
                    .map(|chrom| {
                        println!("{}", chrom);
                        chrom
                    })
                    .collect(),
            )
        } else {
            None
        }
    };
    let (overall_histogram, chrom_to_histogram) = generate_track_histograms(
        &track_filepath,
        num_histogram_buckets.unwrap_or(DEFAULT_HISTOGRAM_NUM_BUCKETS),
        min,
        max,
        bin_size,
        binarize_score,
        filter_chroms,
    )
    .unwrap_or_exit(Some("failed to generate the histogram"));
    let mut keys: Vec<String> =
        chrom_to_histogram.keys().map(|s| s.to_string()).collect();
    keys.sort();
    println!("Per-chromosome histograms:");
    for chrom in keys.iter() {
        println!("Chromosome {}", chrom);
        println!("{}", chrom_to_histogram.get(chrom).unwrap());
    }
    println!("overall_histogram:\n{}", overall_histogram);
}
