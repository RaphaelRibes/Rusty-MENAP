//! Diagnostic: print the spectral probe at a fixed set of thresholds.
//! `cargo run --release --example probe -- 0.80 0.86 0.90`

use std::path::Path;

use rusty_menap::{build_network, load_otu_table, similarity, spectral_test, standardize};

fn main() {
    let thresholds: Vec<f64> = std::env::args()
        .skip(1)
        .filter_map(|a| a.parse().ok())
        .collect();
    let thresholds = if thresholds.is_empty() {
        vec![0.80, 0.82, 0.84, 0.85, 0.86, 0.87, 0.88, 0.90]
    } else {
        thresholds
    };

    let ds = load_otu_table(Path::new("data/Control.txt")).expect("load");
    let (z, _) = standardize(&ds.values).expect("standardize");
    let sim = similarity(&z);

    for t in thresholds {
        let net = build_network(&sim, t);
        match spectral_test(&net.adj, t) {
            Ok(p) => println!(
                "t={:.2} nodes={} edges={} chi2={:.3} df={} p={:.6}",
                t,
                net.len(),
                net.edges(),
                p.chi2,
                p.df,
                p.p
            ),
            Err(e) => println!("t={:.2} nodes={} error: {e}", t, net.len()),
        }
    }
}
