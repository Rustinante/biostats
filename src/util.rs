use std::collections::HashSet;

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
