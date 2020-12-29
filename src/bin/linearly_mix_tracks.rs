use biostats::{
    linear_track_mixture::LinearTrackMixture,
    util::{get_default_human_chrom_inclusion_set, get_weighted_track_paths},
};
use clap::{clap_app, Arg};
use program_flow::{
    argparse::{
        extract_boolean_flag, extract_optional_str_arg, extract_str_arg,
    },
    debug_eprint_named_vars, eprint_named_vars, OrExit,
};

const ONE_BIN_SIZE_STR: &str = "1";

fn main() {
    let mut app = clap_app!(linearly_mix_tracks =>
        (about: "Generates a track in BED format that is the linear \
        combinations of multiple tracks stored in BED format.")
    );
    app = app
        .arg(
            Arg::with_name("weighted_tracks_filepath")
                .long("weighted-tracks")
                .short("f")
                .required(true)
                .takes_value(true)
                .long_help(
                    "Path to a file where each line has two fields, where the \
                    second field is the path to a track in BED format, and the \
                    first field is the weight associated with that track. \
                    For example, a file consisting of two lines:\n\
                    0.5 /path/a.bed\n\
                    0.5 /path/b.bed\n\
                    will produce a track that is the average of a.bed and \
                    b.bed.",
                ),
        )
        .arg(
            Arg::with_name("out_path")
                .long("out-path")
                .short("o")
                .takes_value(true)
                .required(true)
                .help("output file path."),
        )
        .arg(
            Arg::with_name("bin_size")
                .long("bin")
                .takes_value(true)
                .long_help(
                    "Group the base pairs into consecutive bins of size \
                    bin_size, aligned at index 0. For each track, the value \
                    for each bin will be the average of the values for each \
                    base pair in the bin.",
                ),
        )
        .arg(Arg::with_name("binarize_score").long("binarize").help(
            "Each line in the original BED files will contribute a \
                    unit score for the corresponding interval",
        ))
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
            Arg::with_name("default_human_chrom")
                .long("default-human-chrom")
                .short("d")
                .conflicts_with("filter_chrom")
                .help(
                    "Only process chromosomes \
                    chr1, chr2, ... chr22, chrX, chrY.",
                ),
        );
    let matches = app.get_matches();
    let weighted_tracks_filepath =
        extract_str_arg(&matches, "weighted_tracks_filepath");

    let out_path = extract_str_arg(&matches, "out_path");

    let bin_size = extract_optional_str_arg(&matches, "bin_size")
        .unwrap_or_else(|| ONE_BIN_SIZE_STR.to_string())
        .parse::<i64>()
        .unwrap_or_exit(Some(format_args!("failed to parse bin_size")));

    let binarize_score = extract_boolean_flag(&matches, "binarize_score");

    let exclude = extract_optional_str_arg(&matches, "exclude");
    let default_human_chrom =
        extract_boolean_flag(&matches, "default_human_chrom");

    eprint_named_vars!(
        weighted_tracks_filepath,
        out_path,
        bin_size,
        binarize_score,
        default_human_chrom
    );
    debug_eprint_named_vars!(exclude);

    let target_chroms = if default_human_chrom {
        Some(get_default_human_chrom_inclusion_set())
    } else {
        None
    };

    let weighted_bed_files: Vec<(f64, String)> =
        get_weighted_track_paths(&weighted_tracks_filepath)
            .unwrap_or_exit(Some("failed to read the weighted tracks file."));

    let mixture = LinearTrackMixture::create(
        weighted_bed_files,
        bin_size,
        binarize_score,
        exclude,
        target_chroms,
    )
    .unwrap_or_exit(Some("failed to linearly mix the tracks"));

    mixture
        .write_to_bed_file(&out_path)
        .unwrap_or_exit(Some("failed to write to the output file"));
}
