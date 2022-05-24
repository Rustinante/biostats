use crate::{top_k::get_top_k_bin_map, util::get_common_refined_binned_iter};
use math::{
    iter::UnionZip, partition::integer_interval_map::IntegerIntervalMap,
};
use std::collections::HashMap;

pub fn get_top_k_fraction_overlap_ratio(
    map1: &IntegerIntervalMap<f64>,
    map2: &IntegerIntervalMap<f64>,
    top_k_fraction: f64,
    bin_size: i64,
) -> Result<f64, String> {
    let k = get_k(map1, map2, top_k_fraction, bin_size);
    eprintln!(
        "=> top {} fraction corresponds to {} bins",
        top_k_fraction, k
    );
    get_top_k_overlap_ratio(map1, map2, k, bin_size)
}

fn get_k(
    map1: &IntegerIntervalMap<f64>,
    map2: &IntegerIntervalMap<f64>,
    top_k_fraction: f64,
    bin_size: i64,
) -> i64 {
    let count =
        get_common_refined_binned_iter(&map1, &map2, bin_size).count() as f64;
    (count * top_k_fraction) as i64
}

pub fn get_top_k_overlap_ratio(
    map1: &IntegerIntervalMap<f64>,
    map2: &IntegerIntervalMap<f64>,
    k: i64,
    bin_size: i64,
) -> Result<f64, String> {
    let top_k_1 = get_top_k_bin_map(map1, k, bin_size)?;
    let top_k_2 = get_top_k_bin_map(map2, k, bin_size)?;

    let iter = get_common_refined_binned_iter(&top_k_1, &top_k_2, bin_size);

    let mut count = 0i64;
    let mut num_overlapped_bins = 0i64;
    for (_interval, values) in iter {
        count += 1;
        if values[0].is_some() && values[1].is_some() {
            num_overlapped_bins += 1;
        }
    }

    Ok((num_overlapped_bins as f64) / (count as f64))
}

pub fn get_top_k_fraction_overlap_ratio_across_chroms(
    chrom_to_int_interval_map_1: &HashMap<String, IntegerIntervalMap<f64>>,
    chrom_to_int_interval_map_2: &HashMap<String, IntegerIntervalMap<f64>>,
    top_k_fraction: f64,
    bin_size: i64,
) -> Result<f64, String> {
    let mut total_num_bins = 0i64;
    let mut num_overlapped_bins = 0i64;

    let union_zipped_iter = chrom_to_int_interval_map_1
        .union_zip(&chrom_to_int_interval_map_2)
        .into_iter();

    let empty_interval_map = IntegerIntervalMap::new();
    for (_chrom, map_list) in union_zipped_iter {
        let map1 = map_list[0].unwrap_or_else(|| &empty_interval_map);
        let map2 = map_list[1].unwrap_or_else(|| &empty_interval_map);
        let k = get_k(&map1, &map2, top_k_fraction, bin_size);
        let top_k_1 = get_top_k_bin_map(map1, k, bin_size)?;
        let top_k_2 = get_top_k_bin_map(map2, k, bin_size)?;
        let iter = get_common_refined_binned_iter(&top_k_1, &top_k_2, bin_size);

        for (_interval, values) in iter {
            total_num_bins += 1;
            if values[0].is_some() && values[1].is_some() {
                num_overlapped_bins += 1;
            }
        }
    }

    Ok((num_overlapped_bins as f64) / (total_num_bins as f64))
}

#[cfg(test)]
mod tests {
    use crate::{
        check_chrom, test_util::create_temp_bed,
        top_k_overlap::get_top_k_overlap_ratio, util::get_chrom_interval_map,
    };
    use biofile::{bed::Bed, util::TrackVariant};
    use math::iter::UnionZip;
    use std::{
        collections::HashMap,
        io::{BufWriter, Write},
    };
    use tempfile::NamedTempFile;

    #[test]
    fn test_top_k_bed() {
        let bed_1_path = create_temp_bed(
            "chr1 100 200 name_1 10\n\
            chr1 200 250 name_2 75\n\
            chr1 300 350 name_2 125\n\
            chr1 400 450 name_2 25\n\
            chr1 450 550 name_2 500\n\
            chr3 2000 2100 name_6 25\n",
        )
        .unwrap();

        let bed_2_path = create_temp_bed(
            "chr1 100 150 name_1 10\n\
            chr1 300 350 name_2 125\n\
            chr1 400 450 name_2 25\n\
            chr1 600 650 name_2 25\n\
            chr3 2000 2100 name_6 25\n",
        )
        .unwrap();

        let chrom_to_interval_map_1 = get_chrom_interval_map(
            &TrackVariant::Bed(Bed::new(bed_1_path.to_str().unwrap(), false)),
            None,
        )
        .unwrap();

        let chrom_to_interval_map_2 = get_chrom_interval_map(
            &TrackVariant::Bed(Bed::new(bed_2_path.to_str().unwrap(), false)),
            None,
        )
        .unwrap();

        let k = 20i64;
        let bin_size = 25;

        let chrom_to_overlap_ratio: HashMap<String, f64> =
            chrom_to_interval_map_1
                .union_zip(&chrom_to_interval_map_2)
                .into_iter()
                .map(|(chrom, map_list)| {
                    (
                        chrom.to_string(),
                        get_top_k_overlap_ratio(
                            &map_list[0].unwrap(),
                            &map_list[1].unwrap(),
                            k,
                            bin_size,
                        )
                        .unwrap(),
                    )
                })
                .collect();

        assert_eq!(chrom_to_overlap_ratio["chr1"], 150f64 / 400f64);
        assert_eq!(chrom_to_overlap_ratio["chr3"], 1f64);
    }
}
