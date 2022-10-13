use clap::Parser;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use std::collections::BinaryHeap;

use std::cmp::Ordering;
use std::cmp::PartialOrd;

use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;

use rayon::prelude::*;

// use std::collections::HashMap;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    #[clap(parse(from_os_str))]
    query: PathBuf,

    #[clap(parse(from_os_str))]
    matchlist: PathBuf,
}

fn check_compatible_downsample(
    me: &KmerMinHash,
    other: &KmerMinHash,
) -> Result<(), sourmash::Error> {
    /*
    if self.num != other.num {
        return Err(Error::MismatchNum {
            n1: self.num,
            n2: other.num,
        }
        .into());
    }
    */
    use sourmash::Error;

    if me.ksize() != other.ksize() {
        return Err(Error::MismatchKSizes);
    }
    if me.hash_function() != other.hash_function() {
        // TODO: fix this error
        return Err(Error::MismatchDNAProt);
    }
    if me.max_hash() < other.max_hash() {
        return Err(Error::MismatchScaled);
    }
    if me.seed() != other.seed() {
        return Err(Error::MismatchSeed);
    }
    Ok(())
}

fn prepare_query(search_sig: &Signature, template: &Sketch) -> Option<KmerMinHash> {
    let mut search_mh = None;
    if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
        search_mh = Some(mh.clone());
    } else {
        // try to find one that can be downsampled
        if let Sketch::MinHash(template_mh) = template {
            for sketch in search_sig.sketches() {
                if let Sketch::MinHash(ref_mh) = sketch {
                    if check_compatible_downsample(&ref_mh, template_mh).is_ok() {
                        let max_hash = max_hash_for_scaled(template_mh.scaled());
                        let mh = ref_mh.downsample_max_hash(max_hash).unwrap();
                        return Some(mh);
                    }
                }
            }
        }
    }
    search_mh
}

struct PrefetchResult {
    name: String,
    minhash: KmerMinHash,
    containment: u64,
}

impl Ord for PrefetchResult {
    fn cmp(&self, other: &PrefetchResult) -> Ordering {
        self.containment.cmp(&other.containment)
    }
}

impl PartialOrd for PrefetchResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for PrefetchResult {
    fn eq(&self, other: &Self) -> bool {
        self.containment == other.containment
    }
}

impl Eq for PrefetchResult {}

fn prefetch(
    query: &KmerMinHash,
    sketchlist: BinaryHeap<PrefetchResult>,
) -> BinaryHeap<PrefetchResult> {
    sketchlist
        .into_par_iter()
        .filter_map(|result| {
            let mut mm = None;
            let searchsig = &result.minhash;
            let containment = searchsig.count_common(query, false);
            if let Ok(containment) = containment {
                if containment > 0 {
                    let result = PrefetchResult {
                        containment,
                        ..result
                    };
                    mm = Some(result);
                }
            }
            mm
        })
        .collect()
}

fn do_countergather<P: AsRef<Path> + std::fmt::Debug>(
    query_filename: P,
    matchlist: P,
) -> Result<(), Box<dyn std::error::Error>> {
    let max_hash = max_hash_for_scaled(100000_u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(31_u32)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    println!("Loading query");
    let mut query = {
        let sigs = Signature::from_path(dbg!(query_filename)).unwrap();

        let mut mm = None;
        for sig in &sigs {
            if let Some(mh) = prepare_query(sig, &template) {
                mm = Some(mh.clone());
                // doesn't this pick the last one to match the template:
                // hmm. @CTB
            }
        }
        mm
    }
    .unwrap();

    println!("Loading matchlist");
    let matchlist_file = BufReader::new(File::open(matchlist)?);

    // build the list of paths to match against.
    let matchlist_paths: Vec<PathBuf> = matchlist_file
        .lines()
        .filter_map(|line| {
            let line = line.unwrap();
            if !line.is_empty() {
                // skip empty lines
                let mut path = PathBuf::new();
                path.push(line);
                Some(path)
            } else {
                None
            }
        })
        .collect();

    // load the sketches in parallel; keep only those with some match.
    let matchlist: BinaryHeap<PrefetchResult> = matchlist_paths
        .par_iter()
        .filter_map(|m| {
            let sigs = Signature::from_path(m).unwrap();

            let mut mm = None;
            for sig in &sigs {
                if let Some(mh) = prepare_query(sig, &template) {
                    if let Ok(containment) = mh.count_common(&query, false) {
                        if containment > 0 {
                            let result = PrefetchResult {
                                name: sig.name(),
                                minhash: mh,
                                containment,
                            };
                            mm = Some(result);
                            break;
                        }
                    }
                }
            }
            mm
        })
        .collect();

    if matchlist.is_empty() {
        println!("No matchlist signatures loaded, exiting.");
        return Ok(());
    }

    let mut matching_sketches = matchlist;

    // loop until no more matching sketches -
    while !matching_sketches.is_empty() {
        println!("remaining: {} {}", query.size(), matching_sketches.len());
        let best_element = matching_sketches.peek().unwrap();

        // remove!
        println!("removing {}", best_element.name);
        query.remove_from(&best_element.minhash)?;

        // recalculate remaining containments between query and all sketches.
        matching_sketches = prefetch(&query, matching_sketches);
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let opts = Cli::parse();

    do_countergather(opts.query, opts.matchlist)
}
