use biostats::{refined_bed_zipper::RefinedBedZipper, util::get_track_paths};
use clap::{clap_app, Arg};
use program_flow::{
    argparse::{
        extract_numeric_arg, extract_optional_numeric_arg, extract_str_arg,
    },
    eprint_named_vars, OrExit,
};

fn main() {
    let mut app = clap_app!(zip_refined_beds =>
        (about: "Given N binned outputs of refine_bed with the same bin size, \
        zips the values into a single file where each line is of the form \
        chrom start end_exclusive value_1 value_2 ... value_N")
    );
    app = app
        .arg(
            Arg::with_name("refined_beds_path_file")
                .takes_value(true)
                .required(true)
                .help(
                    "Path to a file in which each line is the path \
                    to a refined BED file",
                ),
        )
        .arg(
            Arg::with_name("out_path")
                .takes_value(true)
                .required(true)
                .help("output file path"),
        )
        .arg(
            Arg::with_name("interval_length")
                .long("interval-length")
                .short("l")
                .takes_value(true)
                .required(true)
                .help("Each interval in each BED file must have this size"),
        )
        .arg(
            Arg::with_name("alignment")
                .long("alignment")
                .short("a")
                .takes_value(true)
                .long_help(
                    "The start coordinate of each interval must be a multiple \
                    of the alignment, default to 0.",
                ),
        )
        .arg(
            Arg::with_name("default_value")
                .long("default-value")
                .short("d")
                .help(
                    "If a BED is missing the value corresponding to an \
                    interval, this default value will be used as the value \
                    for that interval for the BED file.",
                ),
        );
    let matches = app.get_matches();
    let refined_beds_path_file =
        extract_str_arg(&matches, "refined_beds_path_file");

    let out_path = extract_str_arg(&matches, "out_path");

    let interval_length: i64 = extract_numeric_arg(&matches, "interval_length")
        .unwrap_or_exit(Some(format!("failed to parse --interval-length")));

    let alignment: i64 = extract_optional_numeric_arg(&matches, "alignment")
        .unwrap_or_exit(Some(format!("failed to parse --alignment")))
        .unwrap_or(0);

    let default_value: f64 =
        extract_optional_numeric_arg(&matches, "default_value")
            .unwrap_or_exit(Some(format!("failed to parse --default-value")))
            .unwrap_or(0.);

    let refined_bed_paths: Vec<String> =
        get_track_paths(&refined_beds_path_file)
            .unwrap_or_exit(Some("failed to get the list of BED file paths"));

    eprint_named_vars!(
        refined_beds_path_file,
        out_path,
        interval_length,
        alignment,
        default_value
    );

    let zipper = RefinedBedZipper::new(
        refined_bed_paths,
        alignment,
        interval_length,
        default_value,
    );
    zipper
        .write_to_file(&out_path)
        .unwrap_or_exit(Some(format!("failed to write to {}", out_path)));
}
