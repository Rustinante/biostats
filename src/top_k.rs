use math::{
    interval::I64Interval,
    iter::{AggregateOp, IntoBinnedIntervalIter},
    partition::integer_interval_map::IntegerIntervalMap,
};
use std::{
    cmp::{Ord, Ordering, PartialOrd, Reverse},
    collections::BinaryHeap,
};

#[derive(PartialEq)]
struct HeapItem {
    interval: I64Interval,
    val: f64,
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.val.partial_cmp(&other.val)
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(&other).expect("failed to compare float")
    }
}

impl Eq for HeapItem {}

pub fn get_top_k(
    interval_map: &IntegerIntervalMap<f64>,
    k: i64,
    bin_size: i64,
) -> Result<IntegerIntervalMap<f64>, String> {
    let binned_iter = interval_map.iter().into_binned_interval_iter(
        bin_size,
        AggregateOp::Average,
        Box::new(|item| (*item.0, *item.1)),
    );

    let mut heap = BinaryHeap::new();

    for (interval, val) in binned_iter {
        heap.push(Reverse(HeapItem {
            interval,
            val,
        }));
        if (heap.len() as i64) > k {
            heap.pop();
        }
    }
    let mut top_interval_map = IntegerIntervalMap::new();
    for item in heap.into_iter() {
        top_interval_map.aggregate(item.0.interval, item.0.val);
    }
    Ok(top_interval_map)
}

#[cfg(test)]
mod tests {
    use crate::{check_chrom, top_k::get_top_k, util::get_chrom_interval_map};
    use biofile::{bed::Bed, util::TrackVariant};
    use math::{
        interval::I64Interval,
        partition::integer_interval_map::IntegerIntervalMap,
    };
    use std::{
        collections::HashMap,
        io::{BufWriter, Write},
    };
    use tempfile::NamedTempFile;

    #[test]
    fn test_top_k_bed() {
        let bed_1_path = {
            let bed_1 = NamedTempFile::new().unwrap();
            {
                let mut writer = BufWriter::new(&bed_1);
                writer
                    .write_fmt(format_args!(
                        "chr1 100 200 name_1 10\n\
                        chr1 200 250 name_2 75\n\
                        chr1 300 350 name_2 125\n\
                        chr1 400 450 name_2 25\n\
                        chr1 450 500 name_2 500\n\
                        chr3 2000 2100 name_6 25\n"
                    ))
                    .unwrap();
            }
            bed_1.into_temp_path()
        };
        let chrom_to_interval_map = get_chrom_interval_map(
            &TrackVariant::Bed(Bed::new(bed_1_path.to_str().unwrap(), false)),
            None,
        )
        .unwrap();

        let top_k_map: HashMap<String, IntegerIntervalMap<f64>> =
            chrom_to_interval_map
                .iter()
                .map(|(chrom, interval_map)| {
                    let top_k = get_top_k(&interval_map, 3, 50).unwrap();
                    (chrom.to_string(), top_k)
                })
                .collect();

        assert_eq!(top_k_map.len(), 2);
        assert_eq!(top_k_map["chr1"].len(), 3);
        assert_eq!(top_k_map["chr3"].len(), 2);

        {
            let mut chr1_map_iter = top_k_map["chr1"].iter();
            check_chrom!(
                chr1_map_iter,
                (200, 249, 75.),
                (300, 349, 125.),
                (450, 499, 500.)
            );
        }

        {
            let mut chr3_map_iter = top_k_map["chr3"].iter();
            check_chrom!(chr3_map_iter, (2000, 2049, 25.), (2050, 2099, 25.));
        }
    }
}
