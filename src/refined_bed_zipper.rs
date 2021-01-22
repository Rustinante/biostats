use biofile::bed::{Bed, BedDataLine, BedDataLineIter};
use math::traits::ToIterator;
use std::{
    collections::HashSet,
    fs::OpenOptions,
    io::{BufWriter, Write},
};

type Chrom = String;
type Coord = i64;
type Value = f64;

#[derive(Clone, Debug, PartialEq)]
pub struct RefinedBedZipper {
    // the output from refine_bed where the intervals are binned and the
    // chromosomes and intervals are sorted
    refined_bed_paths: Vec<String>,

    // the start coordinate of each interval must be divisible by this
    // alignment.
    alignment: Coord,

    // end_exclusive - start for each line in the BED file
    interval_length: Coord,

    default_value: Value,
}

impl RefinedBedZipper {
    pub fn new(
        refined_bed_paths: Vec<String>,
        alignment: Coord,
        interval_length: Coord,
        default_value: Value,
    ) -> RefinedBedZipper {
        RefinedBedZipper {
            refined_bed_paths,
            alignment,
            interval_length,
            default_value,
        }
    }

    pub fn write_to_file(
        &self,
        out_path: &str,
    ) -> Result<(), biofile::error::Error> {
        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(out_path)?;
        let mut writer = BufWriter::new(file);

        for ZippedBedGraphLine {
            chrom,
            start,
            end_exclusive,
            values,
        } in self.try_to_iter()?
        {
            // note that the end coordinate is exclusive in the BED format
            write!(&mut writer, "{}\t{}\t{}", chrom, start, end_exclusive,)?;
            for v in values.iter() {
                write!(&mut writer, "\t{}", v)?;
            }
            writeln!(&mut writer,)?;
        }
        Ok(())
    }
}

pub trait TryToIter<'s, I: Iterator<Item = R>, R, E> {
    fn try_to_iter(&'s self) -> Result<I, E>;
}

impl<'s> TryToIter<'s, RefinedBedZipperIter, ZippedBedGraphLine, String>
    for RefinedBedZipper
{
    fn try_to_iter(&'s self) -> Result<RefinedBedZipperIter, String> {
        let bed_reserves: Vec<BedReserve> = self
            .refined_bed_paths
            .iter()
            .map(|p| BedReserve::new(Bed::new(p, false).to_iter()))
            .collect();

        RefinedBedZipperIter::new(
            bed_reserves,
            self.alignment,
            self.interval_length,
            self.default_value,
        )
    }
}

pub struct ZippedBedGraphLine {
    pub chrom: Chrom,
    pub start: Coord,
    pub end_exclusive: Coord,
    pub values: Vec<Value>,
}

struct RefinedBedZipperIter {
    bed_reserves: Vec<BedReserve>,

    // the start coordinate of each interval must be divisible by this
    // alignment.
    alignment: Coord,

    // end_exclusive - start for each line in the BED file
    interval_length: Coord,

    default_value: Value,
}

impl RefinedBedZipperIter {
    fn new(
        bed_reserves: Vec<BedReserve>,
        alignment: Coord,
        interval_length: Coord,
        default_value: Value,
    ) -> Result<RefinedBedZipperIter, String> {
        if alignment < 0 {
            Err(format!(
                "alignment cannot be negative, received {}",
                alignment
            ))
        } else if interval_length <= 0 {
            Err(format!(
                "interval_legnth must be positive, received {}",
                interval_length
            ))
        } else {
            Ok(RefinedBedZipperIter {
                bed_reserves,
                alignment,
                interval_length,
                default_value,
            })
        }
    }
}

impl Iterator for RefinedBedZipperIter {
    type Item = ZippedBedGraphLine;

    fn next(&mut self) -> Option<Self::Item> {
        let all_chroms: HashSet<Chrom> = self
            .bed_reserves
            .iter()
            .filter_map(|reserve| {
                reserve.current().map(|(chrom, ..)| chrom.to_string())
            })
            .collect();

        let min_chrom = all_chroms.into_iter().min();

        // no minimum chrom means all the bed reserves are exhausted
        match min_chrom {
            None => None,
            Some(min_chrom) => {
                let all_starts: HashSet<Coord> = self
                    .bed_reserves
                    .iter()
                    .filter_map(|reserve| match reserve.current() {
                        None => None,
                        Some((chrom, start, ..)) => {
                            if chrom == &min_chrom {
                                Some(*start)
                            } else {
                                None
                            }
                        }
                    })
                    .collect();

                let min_start = all_starts.into_iter().min().expect(
                    "failed to find the minimum start for the current \
                    chromosome in the bed reserves",
                );

                let alignment = self.alignment;
                let interval_length = self.interval_length;
                let default_value = self.default_value;

                let values: Vec<Value> = self
                    .bed_reserves
                    .iter_mut()
                    .map(|reserve| match reserve.current() {
                        None => default_value,
                        Some((chrom, start, _end_exclusive, value)) => {
                            if chrom == &min_chrom && start == &min_start {
                                let value = if let Some(value) = value {
                                    *value
                                } else {
                                    default_value
                                };
                                reserve
                                    .advance(alignment, interval_length)
                                    .expect("failed to expect the BedReserve");
                                value
                            } else {
                                default_value
                            }
                        }
                    })
                    .collect();

                Some(ZippedBedGraphLine {
                    chrom: min_chrom,
                    start: min_start,
                    end_exclusive: min_start + self.interval_length,
                    values,
                })
            }
        }
    }
}

struct BedReserve {
    bed_iter: BedDataLineIter<Value>,

    // (chrom, start, end_exclusive, value)
    current_chrom_coordinates: Option<(Chrom, Coord, Coord, Option<Value>)>,

    // does not include the current_chrom
    past_chroms: HashSet<Chrom>,
}

impl BedReserve {
    fn new(mut bed_iter: BedDataLineIter<Value>) -> BedReserve {
        let current_chrom_coordinates = match bed_iter.next() {
            None => None,
            Some(BedDataLine {
                chrom,
                start,
                end,
                score,
                ..
            }) => Some((chrom, start, end, score)),
        };

        BedReserve {
            bed_iter,
            current_chrom_coordinates,
            past_chroms: HashSet::new(),
        }
    }

    fn current(&self) -> Option<&(Chrom, Coord, Coord, Option<Value>)> {
        self.current_chrom_coordinates.as_ref()
    }

    fn advance(
        &mut self,
        alignment: Coord,
        interval_length: Coord,
    ) -> Result<Option<&(Chrom, Coord, Coord, Option<Value>)>, String> {
        if let Some(BedDataLine {
            chrom,
            start,
            end,
            score,
            ..
        }) = self.bed_iter.next()
        {
            if self.past_chroms.contains(&chrom) {
                return Err(format!(
                    "different chromosomes cannot interleave, encountered \
                    chromosome {} again: ",
                    &chrom
                ));
            }

            if end - start != interval_length {
                return Err(format!(""));
            }

            if start % interval_length != alignment {
                return Err(format!(""));
            }

            let insert_to_past_chroms = match self.current() {
                None => None,
                Some((old_chrom, ..)) => {
                    if old_chrom != &chrom {
                        Some(old_chrom.clone())
                    } else {
                        None
                    }
                }
            };

            if let Some(old_chrom) = insert_to_past_chroms {
                self.past_chroms.insert(old_chrom);
            }

            if let Some((
                old_chrom,
                _old_start,
                old_end_exclusive,
                _old_value,
            )) = self.current()
            {
                if old_chrom == &chrom {
                    if start < *old_end_exclusive {
                        return Err(format!(
                            "the intervals must be sorted in increasing order, \
                            new start {} < old_end_exclusive {}",
                            start, old_end_exclusive
                        ));
                    }
                } else {
                    if old_chrom > &chrom {
                        return Err(format!(
                            "the chromosomes must be sorted in increasing \
                            order, new chrom {} < old chrom {}",
                            chrom, old_chrom
                        ));
                    }
                }
            }

            self.current_chrom_coordinates = Some((chrom, start, end, score));
        } else {
            // exhausted all lines
            // move the current chrom into the the past_chroms
            match &self.current_chrom_coordinates {
                None => {}
                Some((chrom, ..)) => {
                    self.past_chroms.insert(chrom.to_string());
                }
            }
            self.current_chrom_coordinates = None;
        }

        Ok(self.current())
    }
}
