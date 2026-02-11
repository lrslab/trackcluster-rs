use std::path::PathBuf;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Input BED12 / bigGenePred file
    #[arg(short, long)]
    pub input: PathBuf,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let reader = crate::io::bed::read_bed12(&args.input)?;

    let mut transcript_count: u64 = 0;
    let mut exon_count: u64 = 0;

    for result in reader {
        let transcript = result?;
        transcript_count += 1;
        exon_count += transcript.exons.len() as u64;
    }

    println!("records\t{transcript_count}");
    println!("exons\t{exon_count}");

    Ok(())
}
