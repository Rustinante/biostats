use biostats::track_correlation::{compute_track_correlations, ValueTransform};
use std::{collections::HashSet, path::PathBuf};

const EPSILON: f64 = 1e-8;

macro_rules! assert_almost_eq {
    ($a:expr, $b:expr) => {
        assert!(($a - $b).abs() < EPSILON);
    };
    ($a:expr, $b:expr, $epsilon:expr) => {
        assert!(($a - $b).abs() < $epsilon);
    };
}

macro_rules! assert_vec_almost_eq {
    ($a:expr, $b:expr) => {
        for (x, y) in $a.iter().zip($b.iter()) {
            assert_almost_eq!(x, y, EPSILON);
        }
    };
    ($a:expr, $b:expr, $epsilon:expr) => {
        for (x, y) in $a.iter().zip($b.iter()) {
            assert_almost_eq!(x, y, $epsilon);
        }
    };
}

fn get_bed_path(filename: &str) -> PathBuf {
    let mut path = std::env::current_dir().unwrap();
    path.push(filename);
    path
}

#[test]
fn test_with_identical_track() {
    let chroms: HashSet<String> = vec!["chr1".into()].into_iter().collect();
    let (chrom_correlations, overall_correlations) =
        compute_track_correlations(
            get_bed_path("tests/test_1.bed").to_str().unwrap(),
            get_bed_path("tests/test_2.bed").to_str().unwrap(),
            vec![0, 1, 5, 17],
            chroms,
            ValueTransform::Identity,
        )
        .unwrap();

    assert_eq!(chrom_correlations.len(), 1);
    assert_eq!(chrom_correlations.first().unwrap().0, "chr1");
    assert_vec_almost_eq!(chrom_correlations.first().unwrap().1, vec![
        1., 1., 1., 1.
    ]);

    assert_eq!(overall_correlations.len(), 4);
    assert_vec_almost_eq!(overall_correlations, vec![1., 1., 1., 1.]);
}

#[test]
fn test_with_overlapping_intervals() {
    let chroms: HashSet<String> =
        vec!["chr1".into(), "chr2".into()].into_iter().collect();
    let test_3_bed_path = get_bed_path("tests/test_3.bed");
    let test_4_bed_path = get_bed_path("tests/test_4.bed");
}
