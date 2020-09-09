use biofile::bed::Chrom;
use math::{
    iter::{AsUnionZipped, IntoUnionZip},
    partition::integer_interval_map::IntegerIntervalMap,
};
use num::Num;
use std::{
    collections::{HashMap, HashSet},
    fs::OpenOptions,
    io,
    io::{BufRead, BufReader},
};

/// Returns chr1 to chr22 inclusive, together with chrX and chrY.
pub fn get_default_human_chrom_inclusion_set() -> HashSet<String> {
    let mut chrom_list: HashSet<String> = (1..=22)
        .collect::<Vec<usize>>()
        .iter()
        .map(|i| format!("chr{}", i))
        .collect();
    chrom_list.insert("chrX".to_string());
    chrom_list.insert("chrY".to_string());
    chrom_list
}

pub fn extract_chrom_names(filepath: &str) -> Result<Vec<String>, io::Error> {
    let buf_reader =
        BufReader::new(OpenOptions::new().read(true).open(filepath)?);
    Ok(buf_reader
        .lines()
        .filter_map(|line| {
            let chrom = line.unwrap().trim().to_string();
            if chrom.len() > 0 {
                Some(chrom)
            } else {
                None
            }
        })
        .collect())
}

///
/// * `list_of_chrom_interval_maps`: a vector of maps each mapping chromosomes
///   to integer interval maps.
/// * `target_chroms`: chromosomes not in this set will be ignored.
/// * `default_map`: should be a reference to an empty integer interval map.
pub fn get_union_zipped_chrom_interval_maps<'a, D: Copy + Num>(
    list_of_chrom_interval_maps: Vec<&'a HashMap<Chrom, IntegerIntervalMap<D>>>,
    target_chroms: &HashSet<String>,
    default_map: &'a IntegerIntervalMap<D>,
) -> HashMap<Chrom, Vec<&'a IntegerIntervalMap<D>>> {
    list_of_chrom_interval_maps
        .iter()
        .skip(1)
        .fold(
            list_of_chrom_interval_maps
                .first()
                .expect("there must be at least one chrom interval map")
                .as_union_zipped(),
            |zipped, map| zipped.into_union_zip(map),
        )
        .into_iter()
        .filter_map(|(chrom, map_list)| {
            if target_chroms.contains(&chrom) {
                let maps: Vec<&IntegerIntervalMap<D>> = map_list
                    .into_iter()
                    .map(|map| map.unwrap_or_else(|| default_map))
                    .collect();
                Some((chrom, maps))
            } else {
                None
            }
        })
        .collect()
}

pub fn get_sorted_keys<K: Clone + Ord, V>(map: &HashMap<K, V>) -> Vec<K> {
    let mut keys: Vec<K> = map.keys().cloned().collect();
    keys.sort();
    keys
}
