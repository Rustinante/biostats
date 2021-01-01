use biofile::bed::{Bed, BedDataLine, BedDataLineIter, BedWriter, Chrom};
use math::{
    interval::{traits::Interval, I64Interval},
    iter::{AggregateOp, IntoBinnedIntervalIter},
    partition::integer_interval_map::IntegerIntervalMap,
    set::traits::Set,
    traits::ToIterator,
};
use std::collections::{HashMap, HashSet};

pub struct TrackBinnerStats {
    pub num_duplicate_lines: Option<i64>,
}

pub fn write_common_refined_binned_track(
    track_filepath: &str,
    out_path: &str,
    bin_size: i64,
    unique: bool,
    binarize_score: bool,
    filter_chroms: Option<HashSet<String>>,
    debug: bool,
) -> Result<TrackBinnerStats, biofile::error::Error> {
    let mut visited = HashSet::new();
    let mut num_pcr_duplicates = 0i64;

    let mut chrom_to_interval_map =
        HashMap::<Chrom, IntegerIntervalMap<f64>>::new();

    let bed = Bed::new(track_filepath, binarize_score);
    for BedDataLine {
        chrom,
        start,
        end,
        name: _,
        score,
        strand,
    } in bed.to_iter(): BedDataLineIter<f64>
    {
        if filter_chroms.is_some()
            && !filter_chroms.as_ref().unwrap().contains(&chrom)
        {
            continue;
        }
        if unique {
            if !visited.insert((chrom.clone(), start, end, strand)) {
                // duplicate PCR reads
                num_pcr_duplicates += 1;
                if debug {
                    eprintln!(
                        "PCR duplciate (chrom, start, end, strand): \
                            ({}, {}, {}, {:?})",
                        chrom, start, end, strand
                    )
                }
                continue;
            }
        }
        assert!(
            end > 0,
            "the end coordinate must be positive, encountered \
                (chrom, start, end): ({}, {}, {})",
            chrom,
            start,
            end,
        );

        let interval_map = chrom_to_interval_map
            .entry(chrom)
            .or_insert_with(IntegerIntervalMap::new);

        interval_map
            .aggregate(I64Interval::new(start, end - 1), score.unwrap_or(0.));
    }

    let mut writer = BedWriter::new(out_path)?;
    for chrom in crate::util::get_sorted_keys(&chrom_to_interval_map) {
        let mut data_line = BedDataLine {
            chrom: chrom.to_string(),
            start: 0,
            end: 0,
            name: None,
            score: None::<f64>,
            strand: None,
        };
        let interval_map = &chrom_to_interval_map[&chrom];

        let write_bed_line = |(interval, value): (I64Interval, f64)|
            -> Result<(), biofile::error::Error>{
            if !interval.is_empty() {
                data_line.start = interval.get_start();
                // the end is exclusive in the BED format
                data_line.end = interval.get_end() + 1i64;
                data_line.score = Some(value);
                writer.write_bed_line(&data_line)?;
            }
            Ok(())
        };

        if bin_size == 0 {
            interval_map
                .iter()
                .map(|(&interval, &val)| (interval, val))
                .try_for_each(write_bed_line)?;
        } else {
            interval_map
                .iter()
                .into_binned_interval_iter(
                    bin_size,
                    AggregateOp::Average,
                    Box::new(|item| (*item.0, *item.1)),
                )
                .try_for_each(write_bed_line)?;
        }
    }
    Ok(TrackBinnerStats {
        num_duplicate_lines: if unique {
            Some(num_pcr_duplicates)
        } else {
            None
        },
    })
}
