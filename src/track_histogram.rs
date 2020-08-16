use biofile::bed::Bed;
use math::{
    histogram::Histogram,
    iter::{AggregateOp, IntoBinnedIntervalIter},
    traits::Collecting,
};
use std::collections::{HashMap, HashSet};

type Chrom = String;
type Value = f64;

/// `bed_track_filepath`: BED file path.
/// `histogram_num_buckets`: the number of buckets evenly partitioning the range `[min, max]`
///                          for the histogram.
/// `min`: minimum value for the histogram.
/// `max`: maximum value for the histogram.
/// `bin_size`: the bin size with which to bin the genomic track. The histogram will aggergate
///             the average values of each of the bins. Bin coordinates not overlapping the BED
///             file coordinates are not considered.
pub fn generate_track_histograms(
    bed_track_filepath: &str,
    histogram_num_buckets: usize,
    min: f64,
    max: f64,
    bin_size: i64,
    filter_chroms: Option<HashSet<Chrom>>,
) -> Result<(Histogram<Value>, HashMap<Chrom, Histogram<Value>>), String> {
    let bed = Bed::new(bed_track_filepath);
    let chrom_interval_map = match bed.get_chrom_to_interval_to_val::<Value, _>() {
        Ok(map) => map,
        Err(why) => {
            return Err(format!(
                "failed to get chrom interval map for {}: {}",
                bed_track_filepath, why
            ))
        }
    };
    let mut overall_histogram = Histogram::new(None, histogram_num_buckets, min, max).unwrap();
    let bin_size_f64 = bin_size as Value;
    let mut chrom_to_histogram: HashMap<Chrom, Histogram<Value>> = HashMap::new();
    let mut keys: Vec<Chrom> = chrom_interval_map.keys().map(|k| k.to_string()).collect();
    keys.sort();
    for chrom in keys.iter() {
        if let Some(filter) = &filter_chroms {
            if !filter.contains(chrom) {
                continue;
            }
        }
        let mut chrom_histogram = Histogram::new(None, histogram_num_buckets, min, max).unwrap();
        let interval_map = chrom_interval_map.get(chrom).unwrap();
        for (_interval, v) in interval_map.iter().into_binned_interval_iter(
            bin_size,
            AggregateOp::Sum,
            Box::new(|item| (*item.0, *item.1)),
        ) {
            let average = v / bin_size_f64;
            overall_histogram.collect(average);
            chrom_histogram.collect(average);
        }
        chrom_to_histogram.insert(chrom.clone(), chrom_histogram);
    }
    Ok((overall_histogram, chrom_to_histogram))
}
