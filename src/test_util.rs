use std::io::{BufWriter, Write};
use tempfile::{NamedTempFile, TempPath};
#[macro_export]
macro_rules! check_chrom {
    (
        $iter:expr, $($start_end_val:expr), *
    ) => {
        {
            $(
                let start = $start_end_val.0;
                let end = $start_end_val.1;
                let val = $start_end_val.2;
                assert_eq!(
                    $iter.next(),
                    Some((&I64Interval::new(start, end), &val))
                );
            )*
        }
        assert_eq!(
            $iter.next(),
            None
        );
    };
}

pub fn create_temp_bed(content: &str) -> std::io::Result<TempPath> {
    let file = NamedTempFile::new().unwrap();
    {
        let mut writer = BufWriter::new(&file);
        writer.write_fmt(format_args!("{}", content))?;
    }
    Ok(file.into_temp_path())
}
