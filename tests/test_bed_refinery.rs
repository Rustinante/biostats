#![feature(type_ascription)]

use biofile::{
    bed::{Bed, BedDataLineIter},
    bedgraph::{BedGraph, BedGraphDataLineIter},
};
use biostats::{
    assert_almost_eq, bed_refinery::BedRefinery, util::manifest_path_join,
};
use math::traits::ToIterator;
use tempfile::NamedTempFile;

#[test]
fn test_unique() {
    let refinery = BedRefinery::<f64>::new(
        manifest_path_join("tests/test_7.bed").to_str().unwrap(),
        true,
        true,
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
        false,
        None,
        false,
    );
    let out_file = NamedTempFile::new().unwrap();
    let out_temp_path = out_file.into_temp_path();
    let out_path = out_temp_path.to_str().unwrap().to_string();
    refinery
        .write_refined_bed(&out_path, 0, false, None, true)
        .unwrap();

    let bedgraph = BedGraph::new(&out_path);

    let expected = vec![
        (2, 3, 1.),
        (3, 4, 8.),
        (4, 5, 7.),
        (5, 8, 9.),
        (8, 10, 2.),
        (10, 21, 12.),
        (21, 40, 2.),
    ];
    assert_eq!(
        (bedgraph.to_iter(): BedGraphDataLineIter<f32>).count(),
        expected.len()
    );

    for ((start, end_exclusive, value), line) in expected
        .into_iter()
        .zip(bedgraph.to_iter(): BedGraphDataLineIter<f32>)
    {
        assert_eq!(line.start, start);
        assert_eq!(line.end_exclusive, end_exclusive);
        assert_eq!(line.value, value)
    }
}

// TODO: test with different bin sizes
// TODO: test with scaling and normalize
