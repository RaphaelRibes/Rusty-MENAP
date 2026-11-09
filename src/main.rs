//! Command-line entry point for Rusty-MENAP.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use rusty_menap::{
    ALPHA, COARSE_STEP, FINE_STEP, MenapError, Outcome, Probe, SCAN_MAX, SEED, ScanConfig,
    build_network, load_otu_table, scan, similarity, standardize,
};

const VERSION: &str = env!("CARGO_PKG_VERSION");

const HELP: &str = "\
rusty-menap — RMT threshold selection for molecular ecological networks

Implements the threshold-selection procedure of Deng et al. 2012
(BMC Bioinformatics 13:113): standardise OTU abundances, take |Pearson r|,
scan candidate thresholds, and keep the first one whose adjacency spectrum is
statistically indistinguishable from a Poisson (random) one.

USAGE:
    rusty-menap [OPTIONS]

OPTIONS:
    -i, --input <PATH>   OTU table, rows = OTUs, columns = samples [default: data/Control.txt]
    -o, --out <DIR>      output directory [default: out]
        --alpha <F>      significance level of the chi-square GOF [default: 0.01]
        --st <F>         use this fixed threshold instead of scanning
        --min <F>        lower bound of the scan [default: 0.30]
        --max <F>        upper bound of the scan [default: 0.99]
        --coarse <F>     coarse scan step [default: 0.10]
        --fine <F>       fine scan step [default: 0.01]
        --no-curve       do not write the per-threshold scan curve
        --no-write       analyse only, write nothing
    -q, --quiet          print only the final summary
    -h, --help           show this help
    -V, --version        show the version

OUTPUTS (in --out):
    network.tsv      edge list: source, target, weight
    summary.json     parameters, fitted statistics, network size
    scan_curve.tsv   one row per evaluated threshold

EXIT CODES:
    0  success
    1  analysis failed (no threshold, degenerate spectrum, bad input)
    2  usage error
";

struct Args {
    input: PathBuf,
    out: PathBuf,
    alpha: f64,
    fixed: Option<f64>,
    min: f64,
    max: f64,
    coarse: f64,
    fine: f64,
    curve: bool,
    write: bool,
    quiet: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            input: PathBuf::from("data/Control.txt"),
            out: PathBuf::from("out"),
            alpha: ALPHA,
            fixed: None,
            min: SEED,
            max: SCAN_MAX,
            coarse: COARSE_STEP,
            fine: FINE_STEP,
            curve: true,
            write: true,
            quiet: false,
        }
    }
}

enum Parsed {
    Run(Args),
    Help,
    Version,
}

fn parse_args() -> Result<Parsed, String> {
    let mut args = Args::default();
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = |name: &str| -> Result<String, String> {
            it.next().ok_or_else(|| format!("{name} needs a value"))
        };
        match flag.as_str() {
            "-h" | "--help" => return Ok(Parsed::Help),
            "-V" | "--version" => return Ok(Parsed::Version),
            "-q" | "--quiet" => args.quiet = true,
            "--no-curve" => args.curve = false,
            "--no-write" => args.write = false,
            "-i" | "--input" => args.input = PathBuf::from(value("--input")?),
            "-o" | "--out" => args.out = PathBuf::from(value("--out")?),
            "--alpha" => args.alpha = parse_f64("--alpha", &value("--alpha")?)?,
            "--st" => args.fixed = Some(parse_f64("--st", &value("--st")?)?),
            "--min" => args.min = parse_f64("--min", &value("--min")?)?,
            "--max" => args.max = parse_f64("--max", &value("--max")?)?,
            "--coarse" => args.coarse = parse_f64("--coarse", &value("--coarse")?)?,
            "--fine" => args.fine = parse_f64("--fine", &value("--fine")?)?,
            other => return Err(format!("unknown option {other:?} (try --help)")),
        }
    }
    if !(0.0..1.0).contains(&args.alpha) {
        return Err("--alpha must be in (0, 1)".into());
    }
    if let Some(st) = args.fixed
        && !(0.0..=1.0).contains(&st)
    {
        return Err("--st must be in [0, 1]".into());
    }
    if !(0.0..=1.0).contains(&args.min) || !(0.0..=1.0).contains(&args.max) {
        return Err("--min and --max must be in [0, 1]".into());
    }
    if args.max <= args.min {
        return Err("--max must be greater than --min".into());
    }
    Ok(Parsed::Run(args))
}

fn parse_f64(name: &str, raw: &str) -> Result<f64, String> {
    match raw.parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(v),
        _ => Err(format!("{name}: {raw:?} is not a finite number")),
    }
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(Parsed::Help) => {
            print!("{HELP}");
            return ExitCode::SUCCESS;
        }
        Ok(Parsed::Version) => {
            println!("rusty-menap {VERSION}");
            return ExitCode::SUCCESS;
        }
        Ok(Parsed::Run(args)) => args,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::from(2);
        }
    };

    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> Result<(), MenapError> {
    let ds = load_otu_table(&args.input)?;
    let (n, m) = ds.values.shape();
    if !args.quiet {
        println!(
            "input     {} ({n} OTUs x {m} samples)",
            args.input.display()
        );
    }
    for w in &ds.warnings {
        eprintln!("warning: {w}");
    }

    let (z, constant) = standardize(&ds.values)?;
    if !constant.is_empty() && !args.quiet {
        println!(
            "constant  {} OTU(s) with zero variance, no correlations",
            constant.len()
        );
    }
    let sim = similarity(&z);

    let (outcome, mode) = match args.fixed {
        Some(st) => {
            let net = build_network(&sim, st);
            if net.is_empty() {
                return Err(MenapError::EmptyNetwork { threshold: st });
            }
            let probe = Probe {
                threshold: st,
                nodes: net.len(),
                edges: net.edges(),
                chi2: f64::NAN,
                df: 0,
                p: f64::NAN,
                accepted: true,
            };
            let outcome = Outcome {
                threshold: st,
                network: net,
                probe,
                at_lower_bound: false,
                full_sweep: false,
                curve: Vec::new(),
            };
            (outcome, "fixed")
        }
        None => {
            let cfg = ScanConfig {
                min: args.min,
                max: args.max,
                coarse: args.coarse,
                fine: args.fine,
                alpha: args.alpha,
            };
            let out = scan(&sim, &cfg)?;
            if !args.quiet {
                println!(
                    "scan      {:.2}..{:.2} coarse {:.2} fine {:.2}, alpha = {}",
                    cfg.min, cfg.max, cfg.coarse, cfg.fine, cfg.alpha
                );
                if out.full_sweep {
                    println!(
                        "note      no coarse candidate was compatible with Poisson; swept the \
                         whole range at {:.2}",
                        cfg.fine
                    );
                }
                println!(
                    "\n  {:>5}  {:>6}  {:>8}  {:>8}  {:>3}  {:>8}  decision",
                    "t", "nodes", "edges", "chi2", "df", "p"
                );
                print_curve(&out.curve);
            }
            (out, "rmt")
        }
    };

    println!(
        "\nthreshold {:.2}  ({} mode{})",
        outcome.threshold,
        mode,
        if outcome.at_lower_bound {
            ", accepted at the scan floor"
        } else {
            ""
        }
    );
    println!(
        "network   {} nodes / {} edges, mean degree {:.2}, {} isolates dropped",
        outcome.network.len(),
        outcome.probe.edges,
        outcome.network.mean_degree(),
        n - outcome.network.len()
    );
    if mode == "rmt" {
        println!(
            "fit       chi2 = {:.2}, df = {}, p = {:.4} > alpha = {}",
            outcome.probe.chi2, outcome.probe.df, outcome.probe.p, args.alpha
        );
        if outcome.at_lower_bound {
            eprintln!(
                "warning: the test accepted at the lower bound {:.2}; the transition may lie \
                 below the scanned range",
                outcome.threshold
            );
        }
    }

    if args.write {
        write_outputs(args, &ds.ids, &outcome, mode)?;
    }
    Ok(())
}

fn print_curve(curve: &[Probe]) {
    for p in curve {
        println!(
            "  {:>5.2}  {:>6}  {:>8}  {:>8.2}  {:>3}  {:>8.4}  {}",
            p.threshold,
            p.nodes,
            p.edges,
            p.chi2,
            p.df,
            p.p,
            if p.accepted { "accept" } else { "reject" }
        );
    }
}

fn write_outputs(
    args: &Args,
    ids: &[String],
    outcome: &Outcome,
    mode: &str,
) -> Result<(), MenapError> {
    std::fs::create_dir_all(&args.out)?;

    let mut edges = String::new();
    let p = outcome.network.len();
    for i in 0..p {
        for j in (i + 1)..p {
            let w = outcome.network.adj[(i, j)];
            if w > 0.0 {
                let _ = writeln!(
                    edges,
                    "{}\t{}\t{:.6}",
                    ids[outcome.network.nodes[i]], ids[outcome.network.nodes[j]], w
                );
            }
        }
    }
    std::fs::write(args.out.join("network.tsv"), edges)?;

    if args.curve && mode == "rmt" {
        let mut curve = String::from("threshold\tnodes\tedges\tchi2\tdf\tp\tdecision\n");
        for q in &outcome.curve {
            let _ = writeln!(
                curve,
                "{:.2}\t{}\t{}\t{:.6}\t{}\t{:.6}\t{}",
                q.threshold,
                q.nodes,
                q.edges,
                q.chi2,
                q.df,
                q.p,
                if q.accepted { "accept" } else { "reject" }
            );
        }
        std::fs::write(args.out.join("scan_curve.tsv"), curve)?;
    }

    let json = format!(
        "{{\n  \"input\": \"{}\",\n  \"mode\": \"{}\",\n  \"alpha\": {},\n  \
         \"threshold\": {:.6},\n  \"at_lower_bound\": {},\n  \"nodes\": {},\n  \
         \"edges\": {},\n  \"mean_degree\": {:.6},\n  \"chi2\": {},\n  \"df\": {},\n  \
         \"p\": {}\n}}\n",
        json_escape(&args.input.display().to_string()),
        mode,
        args.alpha,
        outcome.threshold,
        outcome.at_lower_bound,
        outcome.network.len(),
        outcome.probe.edges,
        outcome.network.mean_degree(),
        if outcome.probe.chi2.is_finite() {
            format!("{:.6}", outcome.probe.chi2)
        } else {
            "null".to_string()
        },
        outcome.probe.df,
        if outcome.probe.p.is_finite() {
            format!("{:.6}", outcome.probe.p)
        } else {
            "null".to_string()
        }
    );
    std::fs::write(args.out.join("summary.json"), json)?;

    println!("\nwrote {}", listing(&args.out));
    Ok(())
}

fn listing(dir: &Path) -> String {
    format!(
        "{0}/network.tsv, {0}/summary.json, {0}/scan_curve.tsv",
        dir.display()
    )
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}
