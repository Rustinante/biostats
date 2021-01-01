use biostats::util::{
    extract_chrom_names, get_default_human_chrom_inclusion_set,
};
use clap::{clap_app, Arg};
use program_flow::{
    argparse::{
        extract_boolean_flag, extract_optional_numeric_arg,
        extract_optional_str_arg, extract_str_arg,
    },
    debug_eprint_named_vars, eprint_named_vars, OrExit,
};

fn main() {
    let mut app = clap_app!(common_refine_bed =>
        (about: "Takes the common refinement of the intervals under each \
        chromosome and aggregates their corresponding values, while optionally \
        binning the basepairs. The output intervals will be sorted and \
        disjoint. The input file must be in the BED format")
    );
    app = app
        .arg(
            Arg::with_name("track_filepath")
                .takes_value(true)
                .help("BED filepath for the track"),
        )
        .arg(
            Arg::with_name("out_path")
                .takes_value(true)
                .help("output file path"),
        )
        .arg(
            Arg::with_name("bin_size")
                .long("bin")
                .short("b")
                .takes_value(true)
                .long_help(
                    "Group the base pairs into consecutive bins of size \
                    bin_size, aligned at index 0. A bin size of 0 means \
                    not to bin.",
                ),
        )
        .arg(Arg::with_name("binarize_score").long("binarize").long_help(
            "Each line in the original BED files will contribute a \
            unit score for the corresponding interval",
        ))
        .arg(
            Arg::with_name("debug")
                .long("debug")
                .help("outputs PCR duplicate reads into stderr"),
        )
        .arg(
            Arg::with_name("default_human_chrom")
                .long("default-human-chrom")
                .short("d")
                .conflicts_with("filter_chrom")
                .help(
                    "Only computes histograms for chromosome \
                    chr1, chr2, ... chr22, chrX, chrY",
                ),
        )
        .arg(
            Arg::with_name("filter_chrom")
                .long("filter-chrom")
                .short("f")
                .takes_value(true)
                .help(
                    "Path to a file in which each line is the name of a \
                    chromosome. If provided, will only compute histograms for \
                    the specified list of chromosome names. The leading and \
                    trailing whitespaces of each line will be ignored",
                ),
        )
        .arg(
            Arg::with_name("unique")
                .short("u")
                .long("unique")
                .long_help(
                    "multiple lines with the same \
                    (chromosome, start, end, strand) information will only be \
                    counted once. This option can be used only if \
                    --binarize is set to true",
                ),
        );
    let matches = app.get_matches();
    let bin_size: i64 = extract_optional_numeric_arg(&matches, "bin_size")
        .unwrap_or_exit(Some(format_args!("failed to parse --bin")))
        .unwrap_or(0);

    let binarize_score = extract_boolean_flag(&matches, "binarize_score");
    let debug = extract_boolean_flag(&matches, "debug");
    let default_human_chrom =
        extract_boolean_flag(&matches, "default_human_chrom");

    let filter_chrom = extract_optional_str_arg(&matches, "filter_chrom");
    let out_path = extract_str_arg(&matches, "out_path");
    let track_filepath = extract_str_arg(&matches, "track_filepath");
    let unique = extract_boolean_flag(&matches, "unique");
    if unique && !binarize_score {
        eprintln!("--unique can only be set when --binarize is set");
        std::process::exit(1);
    }

    eprint_named_vars!(
        bin_size,
        binarize_score,
        debug,
        default_human_chrom,
        out_path,
        track_filepath,
        unique
    );
    debug_eprint_named_vars!(filter_chrom);

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

    let stats =
        biostats::track_common_refinery::write_common_refined_binned_track(
            &track_filepath,
            &out_path,
            bin_size,
            unique,
            binarize_score,
            filter_chroms,
            debug,
        )
        .unwrap_or_exit(Some("failed to bin track"));

    match stats.num_duplicate_lines {
        Some(num_duplicates) => {
            println!("number of duplicate lines: {}", num_duplicates)
        }
        None => {}
    }
}
