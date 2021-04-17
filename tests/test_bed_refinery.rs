#![feature(type_ascription)]

use biofile::{
    bed::{Bed, BedDataLineIter},
    bedgraph::{BedGraph, BedGraphDataLineIter},
};
use biostats::{
    assert_almost_eq, bed_refinery::BedRefinery, util::manifest_path_join,
};
use math::traits::ToIterator;
use num::Float;
use std::{fmt::Debug, str::FromStr};
use tempfile::NamedTempFile;

#[test]
fn test_unique() {
    let refinery = BedRefinery::<f64>::new(
        manifest_path_join("tests/test_7.bed").to_str().unwrap(),
        true,
        None,
        true,
        None,
        None,
        false,
    );
    let out_file = NamedTempFile::new().unwrap();
    let out_temp_path = out_file.into_temp_path();
    let out_path = out_temp_path.to_str().unwrap().to_string();

    let bin_size = 0;
    refinery
        .write_refined_bed(&out_path, bin_size, false, None, false)
        .unwrap();

    let bed = Bed::new(&out_path, false);
    let expected = vec![(0, 5, 2.), (5, 8, 1.)];
    assert_eq!(
        (bed.to_iter(): BedDataLineIter<f64>).count(),
        expected.len()
    );
    for ((start, end, score), line) in expected
        .into_iter()
        .zip(bed.to_iter(): BedDataLineIter<f64>)
    {
        assert_eq!(line.start, start);
        assert_eq!(line.end, end);
        assert_almost_eq!(line.score.unwrap(), score);
    }
}

#[test]
fn test_bedgraph() {
    let refinery = BedRefinery::<f64>::new(
        manifest_path_join("tests/test_4.bed").to_str().unwrap(),
        false,
        None,
        false,
        None,
        None,
        false,
    );
    let out_file = NamedTempFile::new().unwrap();
    let out_temp_path = out_file.into_temp_path();
    let out_path = out_temp_path.to_str().unwrap().to_string();
    refinery
        .write_refined_bed(&out_path, 0, false, None, true)
        .unwrap();

    let bedgraph = BedGraph::new(&out_path, false);

    let expected = vec![
        ("chr1", 2, 3, 1.),
        ("chr1", 3, 4, 8.),
        ("chr1", 4, 5, 7.),
        ("chr1", 5, 8, 9.),
        ("chr1", 8, 10, 2.),
        ("chr1", 10, 21, 12.),
        ("chr1", 21, 40, 2.),
    ];

    compare_bedgraph_output(&bedgraph, &expected);
}

#[test]
fn test_max_len() {
    let refinery = BedRefinery::<f64>::new(
        manifest_path_join("tests/test_8.bed").to_str().unwrap(),
        false,
        Some(500usize),
        false,
        None,
        None,
        false,
    );
    let out_file = NamedTempFile::new().unwrap();
    let out_temp_path = out_file.into_temp_path();
    let out_path = out_temp_path.to_str().unwrap().to_string();
    refinery
        .write_refined_bed(&out_path, 0, false, None, true)
        .unwrap();

    let bedgraph = BedGraph::new(&out_path, false);

    let expected = vec![
        ("chr1", 0, 500, 1.),
        ("chr1", 500, 1000, 1.),
        ("chr3", 30, 50, 2.),
        ("chr4", 4, 10, 1.),
    ];
    compare_bedgraph_output(&bedgraph, &expected)
}

// `expected`: (chrom, start, end_exclusive, value)
fn compare_bedgraph_output<
    Value: Debug + Float + FromStr<Err = E>,
    E: Debug,
>(
    bedgraph: &BedGraph,
    expected: &Vec<(&str, i64, i64, Value)>,
) {
    assert_eq!(
        (bedgraph.to_iter(): BedGraphDataLineIter<Value>).count(),
        expected.len()
    );
    for ((chrom, start, end_exclusive, value), line) in expected
        .into_iter()
        .zip(bedgraph.to_iter(): BedGraphDataLineIter<Value>)
    {
        assert_eq!(&line.chrom, chrom);
        assert_eq!(line.start, *start);
        assert_eq!(line.end_exclusive, *end_exclusive);
        assert_eq!(line.value, *value);
    }
}

// TODO: test with different bin sizes
// TODO: test with scaling and normalize
// TODO: test exclude
