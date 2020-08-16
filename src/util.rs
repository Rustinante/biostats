use std::{
    collections::HashSet,
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
    let buf_reader = BufReader::new(OpenOptions::new().read(true).open(filepath)?);
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
