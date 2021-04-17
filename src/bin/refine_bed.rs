use biostats::{
    bed_refinery::BedRefinery,
    util::{extract_chrom_names, get_default_human_chrom_inclusion_set},
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
        binning the basepairs. The output intervals will be disjoint and \
        increasing. The input file must be in the BED format.")
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
            Arg::with_name("exclude")
                .long("exclude")
                .short("v")
                .takes_value(true)
                .help(
                    "Path to a BED-like file where only the chromosome, start \
                    and end fields are required. Lines from other BED files \
                    that overlap with any of the coordinates in this 'exclude' \
                    file will be ignored.",
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
            Arg::with_name("max_len")
                .long("max-len")
                .takes_value(true)
                .help(
                    "If provided, will ignore lines in which \
                    end - start > max_len. Note that in the BED format, \
                    end is exclusive, so end - start is exactly the \
                    number of basepairs.",
                ),
        )
        .arg(
            Arg::with_name("normalize")
                .short("n")
                .long("normalize")
                .help(
                    "Normalize by dividing the value associated with each \
                    basepair by the sum of all values across the basepairs.",
                ),
        )
        .arg(
            Arg::with_name("scale")
                .short("s")
                .long("scale")
                .takes_value(true)
                .long_help(
                    "Scales the value associated with each basepair. \
                    If --normalize is set, the values will be normalized \
                    before being scaled.",
                ),
        )
        .arg(
            Arg::with_name("unique")
                .short("u")
                .long("unique")
                .long_help(
                    "Multiple lines with the same \
                    (chromosome, start, end, strand) information will only be \
                    counted once. This option can be used only if \
                    --binarize is set to true",
                ),
        )
        .arg(Arg::with_name("out_bedgraph").long("out-bedgraph").help(
            "Output will be in the Begraph format, i.e., each line \
                    will consist of 4 fields, \
                    (chromosome, start, end_exclusive, value)",
        ));
    let matches = app.get_matches();
    let bin_size: i64 = extract_optional_numeric_arg(&matches, "bin_size")
        .unwrap_or_exit(Some(format_args!("failed to parse --bin")))
        .unwrap_or(0);

    let binarize_score = extract_boolean_flag(&matches, "binarize_score");
    let debug = extract_boolean_flag(&matches, "debug");
    let default_human_chrom =
        extract_boolean_flag(&matches, "default_human_chrom");

    let exclude = extract_optional_str_arg(&matches, "exclude");
    let filter_chrom = extract_optional_str_arg(&matches, "filter_chrom");
    let out_path = extract_str_arg(&matches, "out_path");
    let track_filepath = extract_str_arg(&matches, "track_filepath");
    let max_len: Option<usize> =
        extract_optional_numeric_arg(&matches, "max_len")
            .unwrap_or_exit(Some("failed to parse the --max-len argument"));

    let normalize = extract_boolean_flag(&matches, "normalize");
    let scale: Option<f64> = extract_optional_numeric_arg(&matches, "scale")
        .unwrap_or_exit(None::<String>);

    let unique = extract_boolean_flag(&matches, "unique");
    if unique && !binarize_score {
        eprintln!("--unique can only be set when --binarize is set");
        std::process::exit(1);
    }
    let out_bedgraph = extract_boolean_flag(&matches, "out_bedgraph");

    eprint_named_vars!(
        bin_size,
        binarize_score,
        debug,
        default_human_chrom,
        out_path,
        track_filepath,
        normalize,
        unique,
        out_bedgraph
    );
    debug_eprint_named_vars!(exclude, filter_chrom, max_len, scale);

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

    let refinery = BedRefinery::<f64>::new(
        &track_filepath,
        unique,
        max_len,
        binarize_score,
        filter_chroms,
        exclude,
        debug,
    );

    refinery
        .write_refined_bed(&out_path, bin_size, normalize, scale, out_bedgraph)
        .unwrap_or_exit(Some("failed to bin track"));

    match refinery.stats().num_duplicate_lines {
        Some(num_duplicates) => {
            println!("number of duplicate lines: {}", num_duplicates)
        }
        None => {}
    }
}
