use biofile::{bed::Chrom, iter::ToChromIntervalValueIter, util::TrackVariant};
use math::{
    iter::{AsUnionZipped, IntoUnionZip},
    partition::integer_interval_map::IntegerIntervalMap,
    set::ordered_integer_set::OrderedIntegerSet,
};
use num::Num;
use std::{
    collections::{HashMap, HashSet},
    fs::OpenOptions,
    io,
    io::{BufRead, BufReader},
    path::PathBuf,
};

#[macro_export]
macro_rules! assert_almost_eq {
    ($a:expr, $b:expr) => {
        assert!(($a - $b).abs() < 1e-8);
    };
    ($a:expr, $b:expr, $epsilon:expr) => {
        assert!(($a - $b).abs() < $epsilon);
    };
}

#[macro_export]
macro_rules! assert_vec_almost_eq {
    ($a:expr, $b:expr) => {
        assert_eq!($a.len(), $b.len());
        for (x, y) in $a.iter().zip($b.iter()) {
            assert!((x - y).abs() < 1e-8);
        }
    };
    ($a:expr, $b:expr, $epsilon:expr) => {
        assert_eq!($a.len(), $b.len());
        for (x, y) in $a.iter().zip($b.iter()) {
            assert!((x - y).abs() < $epsilon);
        }
    };
}

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

pub fn get_track_paths(
    track_paths_file: &str,
) -> Result<Vec<String>, std::io::Error> {
    let buf_reader =
        BufReader::new(OpenOptions::new().read(true).open(track_paths_file)?);

    Ok(buf_reader
        .lines()
        .filter_map(|line| {
            let path = line.unwrap().trim().to_string();
            if path.is_empty() {
                None
            } else {
                Some(path)
            }
        })
        .collect())
}

pub fn get_weighted_track_paths(
    filepath: &str,
) -> Result<Vec<(f64, String)>, std::io::Error> {
    let buf_reader =
        BufReader::new(OpenOptions::new().read(true).open(filepath)?);

    Ok(buf_reader
        .lines()
        .filter_map(|line| {
            let tokens: Vec<String> = line
                .unwrap()
                .trim()
                .split_whitespace()
                .map(|s| s.to_string())
                .collect();
            if tokens.is_empty() {
                None
            } else {
                assert_eq!(
                    tokens.len(),
                    2,
                    "Each field in the weighted tracks \
                    file must have exactly two fields"
                );
                let weight = tokens[0]
                    .parse::<f64>()
                    .expect("failed to parse the weight");
                Some((weight, tokens[1].clone()))
            }
        })
        .collect())
}

pub fn get_chrom_interval_map(
    track: &TrackVariant,
    exclude: Option<&HashMap<String, OrderedIntegerSet<i64>>>,
) -> Result<HashMap<String, IntegerIntervalMap<f64>>, String> {
    let result = match track {
        TrackVariant::Bed(bed) => {
            ToChromIntervalValueIter::get_chrom_to_interval_to_val(bed, exclude)
        }
        TrackVariant::BedGraph(bedgraph) => {
            ToChromIntervalValueIter::get_chrom_to_interval_to_val(
                bedgraph, exclude,
            )
        }
    };
    result.map_err(|why| {
        format!("failed to get the first chrom interval map: {}", why)
    })
}

///
/// * `list_of_chrom_interval_maps`: a vector of maps each mapping chromosomes
///   to integer interval maps.
/// * `target_chroms`: chromosomes not in this set will be ignored.
/// * `default_map`: should be a reference to an empty integer interval map.
pub fn get_union_zipped_chrom_interval_maps<'a, D: Copy + Num>(
    list_of_chrom_interval_maps: Vec<&'a HashMap<Chrom, IntegerIntervalMap<D>>>,
    target_chroms: Option<&HashSet<String>>,
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
            if let Some(target_chroms) = target_chroms {
                if target_chroms.contains(&chrom) {
                    let maps: Vec<&IntegerIntervalMap<D>> = map_list
                        .into_iter()
                        .map(|map| map.unwrap_or_else(|| default_map))
                        .collect();
                    Some((chrom, maps))
                } else {
                    None
                }
            } else {
                let maps: Vec<&IntegerIntervalMap<D>> = map_list
                    .into_iter()
                    .map(|map| map.unwrap_or_else(|| default_map))
                    .collect();
                Some((chrom, maps))
            }
        })
        .collect()
}

pub fn get_sorted_keys<K: Clone + Ord, V>(map: &HashMap<K, V>) -> Vec<K> {
    let mut keys: Vec<K> = map.keys().cloned().collect();
    keys.sort();
    keys
}

pub fn manifest_path_join(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push(filename);
    path
}
