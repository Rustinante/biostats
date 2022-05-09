use biofile::{bed::Bed, util::TrackVariant};
use biostats::{
    assert_vec_almost_eq, track_correlation::ValueTransform,
    util::manifest_path_join,
};
use std::collections::HashSet;

#[test]
fn test_identical_tracks() {
    let chroms: HashSet<String> = vec!["chr1".into()].into_iter().collect();
    let first_track = TrackVariant::Bed(Bed::new(
        manifest_path_join("tests/test_1.bed").to_str().unwrap(),
        false,
    ));
    let second_track = TrackVariant::Bed(Bed::new(
        manifest_path_join("tests/test_2.bed").to_str().unwrap(),
        false,
    ));
    let (chrom_correlations, overall_correlations) =
        biostats::track_correlation::compute_track_correlations(
            &first_track,
            &second_track,
            &vec![0, 1, 5, 17],
            Some(chroms),
            ValueTransform::Identity,
            None,
            None,
        )
        .unwrap();

    assert_eq!(chrom_correlations.len(), 1);
    assert_eq!(chrom_correlations.first().unwrap().0, "chr1");
    assert_vec_almost_eq!(chrom_correlations.first().unwrap().1, vec![
        1., 1., 1., 1.
    ]);

    assert_vec_almost_eq!(overall_correlations, vec![1., 1., 1., 1.]);
}

#[test]
fn test_single_chrom() {
    let chroms: HashSet<String> = vec!["chr1".into()].into_iter().collect();

    let first_track = TrackVariant::Bed(Bed::new(
        manifest_path_join("tests/test_3.bed").to_str().unwrap(),
        false,
    ));
    let second_track = TrackVariant::Bed(Bed::new(
        manifest_path_join("tests/test_4.bed").to_str().unwrap(),
        false,
    ));

    let (chrom_correlations, overall_correlations) =
        biostats::track_correlation::compute_track_correlations(
            &first_track,
            &second_track,
            &vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13],
            Some(chroms),
            ValueTransform::Identity,
            None,
            None,
        )
        .unwrap();

    assert_eq!(chrom_correlations.len(), 1);
    assert_eq!(chrom_correlations[0].0, "chr1");

    assert_vec_almost_eq!(chrom_correlations[0].1, vec![
        0.08274757829205674,
        0.08274757829205674,
        0.085717287656436,
        0.13009057672572488,
        0.05688179313971669,
        0.09410511434089433,
        0.3575532028481459,
        0.11257975517362048,
        0.4927537474122497,
        0.40482242419857467,
        0.03131867977157447,
        0.5556478168884006,
        0.8735302221990147
    ]);

    assert_vec_almost_eq!(overall_correlations, vec![
        0.08274757829205674,
        0.08274757829205674,
        0.085717287656436,
        0.13009057672572488,
        0.05688179313971669,
        0.09410511434089433,
        0.3575532028481459,
        0.11257975517362048,
        0.4927537474122497,
        0.40482242419857467,
        0.03131867977157447,
        0.5556478168884006,
        0.8735302221990147
    ]);
}

#[test]
fn test_two_chroms() {
    let chroms: HashSet<String> =
        vec!["chr1".into(), "chr2".into()].into_iter().collect();

    let first_track = TrackVariant::Bed(Bed::new(
        manifest_path_join("tests/test_5.bed").to_str().unwrap(),
        false,
    ));
    let second_track = TrackVariant::Bed(Bed::new(
        manifest_path_join("tests/test_6.bed").to_str().unwrap(),
        false,
    ));

    let (chrom_correlations, overall_correlations) =
        biostats::track_correlation::compute_track_correlations(
            &first_track,
            &second_track,
            &vec![0, 1, 2, 5],
            Some(chroms),
            ValueTransform::Identity,
            None,
            None,
        )
        .unwrap();

    assert_eq!(chrom_correlations.len(), 2);
    assert_eq!(chrom_correlations[0].0, "chr1");
    assert_eq!(chrom_correlations[1].0, "chr2");

    assert_vec_almost_eq!(chrom_correlations[0].1, vec![
        -0.29462546563051467,
        -0.29462546563051467,
        -0.32507713972384644,
        1.
    ]);

    assert_vec_almost_eq!(chrom_correlations[1].1, vec![
        -0.20938823948535507,
        -0.20938823948535507,
        -0.03846124818029947,
        1.
    ]);

    assert_vec_almost_eq!(overall_correlations, vec![
        0.005916134918602554,
        0.005916134918602554,
        0.19356634819114077,
        0.6687843872007803
    ]);
}
