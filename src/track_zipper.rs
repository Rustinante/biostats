use biofile::bed::{Bed, Chrom};
use math::{
    interval::I64Interval,
    iter::{AggregateOp, CommonRefinementZip, IntoBinnedIntervalIter},
    partition::integer_interval_map::IntegerIntervalMap,
};
use std::collections::{HashMap, HashSet};

type Value = f64;

pub struct TrackZipper {
    bed_files: Vec<Bed>,
    list_of_chrom_interval_maps: Vec<HashMap<Chrom, IntegerIntervalMap<Value>>>,
}

impl TrackZipper {
    pub fn new(
        bed_files: Vec<Bed>,
        exclude_track_filepath: Option<&str>,
    ) -> Result<Self, biofile::error::Error> {
        let exclude = if let Some(path) = exclude_track_filepath {
            // binarize_score is irrelevant for getting the intervals
            Some(Bed::new(path, false).get_chrom_to_intervals())
        } else {
            None
        };

        let list_of_chrom_interval_maps: Vec<
            HashMap<Chrom, IntegerIntervalMap<Value>>,
        > = bed_files
            .iter()
            .map(|bed| bed.get_chrom_to_interval_to_val(exclude.as_ref()))
            .collect::<Result<
                Vec<HashMap<Chrom, IntegerIntervalMap<Value>>>,
                biofile::error::Error,
            >>()?;

        Ok(TrackZipper {
            bed_files,
            list_of_chrom_interval_maps,
        })
    }

    pub fn num_tracks(&self) -> usize {
        self.list_of_chrom_interval_maps.len()
    }

    pub fn chrom_to_binned_zipped_values(
        &self,
        target_chroms: Option<&HashSet<Chrom>>,
        bin_size: i64,
    ) -> Result<
        HashMap<Chrom, Vec<(I64Interval, Vec<Option<Value>>)>>,
        biofile::error::Error,
    > {
        let empty_interval_map = IntegerIntervalMap::new();
        let union_zipped_chrom_interval_maps: HashMap<
            Chrom,
            Vec<&IntegerIntervalMap<Value>>,
        > = crate::util::get_union_zipped_chrom_interval_maps(
            self.list_of_chrom_interval_maps.iter().collect(),
            target_chroms,
            &empty_interval_map,
        );
        let chroms =
            crate::util::get_sorted_keys(&union_zipped_chrom_interval_maps);

        Ok(chroms
            .into_iter()
            .map(|chrom| {
                let interval_maps = &union_zipped_chrom_interval_maps[&chrom];

                let binned_values: Vec<(I64Interval, Vec<Option<Value>>)> =
                    interval_maps
                        .iter()
                        .skip(1)
                        .fold(
                            interval_maps
                                .first()
                                .expect("interval maps cannot be empty")
                                .iter()
                                .into_binned_interval_iter(
                                    bin_size,
                                    AggregateOp::Average,
                                    Box::new(|item| (*item.0, *item.1)),
                                )
                                .into_common_refinement_zipped(),
                            |common_refinement, map| {
                                common_refinement.common_refinement_flat_zip(
                                    map.iter().into_binned_interval_iter(
                                        bin_size,
                                        AggregateOp::Average,
                                        Box::new(|item| (*item.0, *item.1)),
                                    ),
                                )
                            },
                        )
                        .collect();
                (chrom, binned_values)
            })
            .collect())
    }
}
