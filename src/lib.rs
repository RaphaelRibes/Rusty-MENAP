//! Rusty-MENAP: RMT-based threshold selection for molecular ecological networks.
//!
//! Implements the threshold-selection procedure of Deng et al. 2012,
//! "Molecular ecological network analyses", BMC Bioinformatics 13:113,
//! Methods steps (a)-(h):
//!
//! (a) standardise OTU abundances per OTU across samples,
//! (b) similarity `S = |R|` (absolute Pearson correlation),
//! (c) seed threshold 0.30,
//! (d) adjacency `A_pxp = s_ij if s_ij >= st else 0`, where `p <= n` because
//!     nodes whose whole row and column fall below the threshold are removed,
//! (e) eigenvalues of `A`,
//! (f) unfold them with a Gaussian-broadened cumulative density (a local
//!     average, which is one of the two unfolding methods the paper allows),
//! (g) nearest-neighbour spacings of the unfolded spectrum,
//! (h) chi-square goodness-of-fit against the Poisson NNSD `exp(-d)`;
//!     the first threshold that does not reject Poisson is returned.
//!
//! Deliberate choices that the paper leaves open are marked `CHOICE` and
//! listed in the README.

use std::fmt;
use std::path::Path;

use nalgebra::DMatrix;

/// Seed threshold from paper step (c).
pub const SEED: f64 = 0.30;
/// Upper end of the scan (the paper does not specify one).
pub const SCAN_MAX: f64 = 0.99;
/// Coarse scan step.
pub const COARSE_STEP: f64 = 0.10;
/// Fine scan step around the coarse transition.
pub const FINE_STEP: f64 = 0.01;
/// Significance level from the paper ("critical value at a significant level
/// of 0.01").
pub const ALPHA: f64 = 0.01;
/// Below this many levels the NNSD is not estimable.
pub const MIN_LEVELS: usize = 20;
/// Eigenvalues closer to zero than this are treated as the null space of a
/// disconnected graph rather than as part of the bulk.
const ZERO_EIG: f64 = 1e-9;
/// Gaussian-broadening width, in units of the mean level spacing.
const SIGMA_LEVELS: f64 = 5.0;

// ---------------------------------------------------------------------------
// errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum MenapError {
    Io(String),
    Empty,
    NoData,
    Ragged {
        line: usize,
        expected: usize,
        found: usize,
    },
    BadNumber {
        line: usize,
        column: usize,
        value: String,
    },
    TooFewSamples(usize),
    EmptyNetwork {
        threshold: f64,
    },
    TooFewEigenvalues {
        threshold: f64,
        levels: usize,
    },
    Degenerate {
        threshold: f64,
        what: String,
    },
    NoThreshold {
        max: f64,
        last_p: f64,
        reason: String,
    },
}

impl fmt::Display for MenapError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::Empty => write!(f, "input file is empty"),
            Self::NoData => write!(f, "input file has a header but no data rows"),
            Self::Ragged {
                line,
                expected,
                found,
            } => write!(
                f,
                "line {line}: expected {expected} tab-separated columns, found {found}"
            ),
            Self::BadNumber {
                line,
                column,
                value,
            } => write!(
                f,
                "line {line}, column {column}: cannot parse {value:?} as a number"
            ),
            Self::TooFewSamples(n) => {
                write!(
                    f,
                    "need at least 2 samples to compute a correlation, found {n}"
                )
            }
            Self::EmptyNetwork { threshold } => write!(
                f,
                "network is empty at threshold {threshold:.2}: no pair of OTUs exceeds it"
            ),
            Self::TooFewEigenvalues { threshold, levels } => write!(
                f,
                "only {levels} non-degenerate eigenvalues at threshold {threshold:.2}, \
                 need at least {MIN_LEVELS}"
            ),
            Self::Degenerate { threshold, what } => {
                write!(f, "degenerate spectrum at threshold {threshold:.2}: {what}")
            }
            Self::NoThreshold {
                max,
                last_p,
                reason,
            } => write!(
                f,
                "no threshold in the scan below {max:.2} is compatible with Poisson \
                 (last p = {last_p:.4}; {reason})"
            ),
        }
    }
}

impl std::error::Error for MenapError {}

impl From<std::io::Error> for MenapError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e.to_string())
    }
}

// ---------------------------------------------------------------------------
// step (a): load and standardise
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct Dataset {
    /// OTU identifiers, one per row.
    pub ids: Vec<String>,
    /// Abundance table, OTUs x samples.
    pub values: DMatrix<f64>,
    /// Non-fatal notes (duplicate ids, constant rows).
    pub warnings: Vec<String>,
}

/// Read a `ID<TAB>sample...` table. Blank lines are skipped; empty cells are
/// read as zero, which is the conventional encoding of an unobserved OTU.
/// Any other deviation from the header width is an error rather than a silent
/// zero-fill, so a shifted table cannot quietly become a different matrix.
pub fn load_otu_table(path: &Path) -> Result<Dataset, MenapError> {
    let text = std::fs::read_to_string(path)?;
    let text = text.strip_prefix('\u{feff}').unwrap_or(&text);

    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty());

    let (_, header) = lines.next().ok_or(MenapError::Empty)?;
    let header: Vec<&str> = header.split('\t').map(str::trim).collect();
    let n_samples = header.len().saturating_sub(1);
    if n_samples < 2 {
        return Err(MenapError::TooFewSamples(n_samples));
    }

    let mut ids = Vec::new();
    let mut values = Vec::new();
    for (idx, line) in lines {
        let line_no = idx + 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != header.len() {
            return Err(MenapError::Ragged {
                line: line_no,
                expected: header.len(),
                found: fields.len(),
            });
        }
        ids.push(fields[0].trim().to_string());
        for (col, raw) in fields[1..].iter().enumerate() {
            let raw = raw.trim();
            let v = if raw.is_empty() {
                0.0
            } else {
                raw.parse::<f64>().map_err(|_| MenapError::BadNumber {
                    line: line_no,
                    column: col + 2,
                    value: raw.to_string(),
                })?
            };
            values.push(v);
        }
    }

    if ids.is_empty() {
        return Err(MenapError::NoData);
    }

    let mut warnings = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for id in &ids {
        if !seen.insert(id.as_str()) {
            warnings.push(format!("duplicate OTU identifier {id:?}"));
        }
    }

    let n = ids.len();
    let matrix = DMatrix::from_row_slice(n, n_samples, &values);
    Ok(Dataset {
        ids,
        values: matrix,
        warnings,
    })
}

/// Row-wise z-score, the paper's "standardized relative abundance".
/// Returns the standardised matrix and the indices of constant OTUs, which
/// carry no correlation information and are reported rather than silently
/// turned into a zero row.
pub fn standardize(values: &DMatrix<f64>) -> Result<(DMatrix<f64>, Vec<usize>), MenapError> {
    let (n, m) = values.shape();
    if m < 2 {
        return Err(MenapError::TooFewSamples(m));
    }
    let mut z = DMatrix::zeros(n, m);
    let mut constant = Vec::new();
    for r in 0..n {
        let mean = values.row(r).iter().sum::<f64>() / m as f64;
        let var = values
            .row(r)
            .iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f64>()
            / (m as f64 - 1.0);
        let sd = var.sqrt();
        if sd <= 1e-12 {
            constant.push(r);
            continue;
        }
        for c in 0..m {
            z[(r, c)] = (values[(r, c)] - mean) / sd;
        }
    }
    Ok((z, constant))
}

// ---------------------------------------------------------------------------
// step (b): similarity
// ---------------------------------------------------------------------------

/// `S = |R|`, with the diagonal pinned to 1 and values clamped into [0, 1]
/// (floating point can overshoot 1 by an ulp for identical rows).
pub fn similarity(z: &DMatrix<f64>) -> DMatrix<f64> {
    let m = z.ncols() as f64;
    let mut s = (z * z.transpose()).scale(1.0 / (m - 1.0));
    s.apply(|x| {
        *x = x.abs().min(1.0);
    });
    for i in 0..s.nrows() {
        s[(i, i)] = 1.0;
    }
    s
}

// ---------------------------------------------------------------------------
// step (d): adjacency
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Network {
    /// Row indices into the original OTU table, ascending.
    pub nodes: Vec<usize>,
    /// Weighted adjacency over `nodes`, zero diagonal.
    pub adj: DMatrix<f64>,
}

impl Network {
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Undirected edge count: nonzero entries strictly above the diagonal.
    pub fn edges(&self) -> usize {
        let p = self.adj.nrows();
        let mut e = 0;
        for i in 0..p {
            for j in (i + 1)..p {
                if self.adj[(i, j)] > 0.0 {
                    e += 1;
                }
            }
        }
        e
    }

    pub fn degree(&self, i: usize) -> usize {
        (0..self.adj.ncols())
            .filter(|&j| j != i && self.adj[(i, j)] > 0.0)
            .count()
    }

    pub fn mean_degree(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        (0..self.len()).map(|i| self.degree(i)).sum::<usize>() as f64 / self.len() as f64
    }
}

/// Threshold the similarity matrix and keep only connected nodes, as in paper
/// step (d). The diagonal is zero: a self-loop is not an edge, and leaving it
/// at 1.0 would make every row "non-zero" so no node could ever be pruned.
pub fn build_network(sim: &DMatrix<f64>, threshold: f64) -> Network {
    let n = sim.nrows();
    let keep: Vec<usize> = (0..n)
        .filter(|&i| (0..n).any(|j| j != i && sim[(i, j)] >= threshold))
        .collect();
    let p = keep.len();
    let mut adj = DMatrix::zeros(p, p);
    for (a, &i) in keep.iter().enumerate() {
        for (b, &j) in keep.iter().enumerate() {
            if a != b && sim[(i, j)] >= threshold {
                adj[(a, b)] = sim[(i, j)];
            }
        }
    }
    Network { nodes: keep, adj }
}

// ---------------------------------------------------------------------------
// steps (e)-(h): spectral test
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Probe {
    pub threshold: f64,
    pub nodes: usize,
    pub edges: usize,
    pub chi2: f64,
    pub df: usize,
    pub p: f64,
    pub accepted: bool,
}

/// Unfold eigenvalues with a Gaussian-broadened cumulative density.
///
/// `N_av(e) = sum_j Phi((e - e_j)/sigma)`. Writing the levels below `e` as
/// `1 - Phi((e_j - e)/sigma)` keeps every term in the decaying half of the
/// kernel, so the sum can be truncated at 8 sigma without touching the result
/// (the omitted mass is below 1e-15 per level).
pub fn unfold(evs: &[f64]) -> Vec<f64> {
    let n = evs.len();
    let range = evs[n - 1] - evs[0];
    let sigma = if range > 0.0 {
        SIGMA_LEVELS * range / n as f64
    } else {
        1.0
    };
    let cut = 8.0 * sigma;
    let inv = std::f64::consts::SQRT_2 * sigma;
    let mut out = Vec::with_capacity(n);
    for (k, &e) in evs.iter().enumerate() {
        // Everything at or below e_k contributes 1 - Phi((e_k - e_j)/sigma);
        // everything above contributes Phi((e_j - e_k)/sigma). Both arguments
        // are positive, so erfc -- not 1 - erf -- carries the tail without
        // cancellation.
        let mut sum = (k + 1) as f64;
        let lo = evs.partition_point(|&x| x < e - cut);
        for &ej in &evs[lo..=k] {
            sum -= 0.5 * erfc((e - ej) / inv);
        }
        let hi = evs.partition_point(|&x| x <= e + cut);
        for &ej in &evs[k + 1..hi] {
            sum += 0.5 * erfc((ej - e) / inv);
        }
        out.push(sum);
    }
    out
}

/// Nearest-neighbour spacings, normalised to mean 1.
pub fn spacings(levels: &[f64]) -> Vec<f64> {
    let mut d: Vec<f64> = levels
        .windows(2)
        .map(|w| w[1] - w[0])
        .filter(|x| *x > 0.0)
        .collect();
    if d.is_empty() {
        return d;
    }
    let mean = d.iter().sum::<f64>() / d.len() as f64;
    for x in &mut d {
        *x /= mean;
    }
    d
}

/// Chi-square goodness-of-fit of spacings against `exp(-d)`.
///
/// Bins are `round(sqrt(n))` wide and equal-width up to the largest observed
/// spacing; the last bin is open-ended. Adjacent bins are merged until each
/// expectation reaches 5, so no cell can dominate the statistic. `df` is
/// `bins - 2` because the mean spacing is estimated from the same data.
fn chi_square_poisson(spacings: &[f64]) -> Result<(f64, usize, f64), String> {
    let n = spacings.len();
    let max_s = spacings.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max_s.is_nan() || max_s <= 0.0 || max_s.is_infinite() {
        return Err("all spacings are zero or non-finite".into());
    }
    let bins = ((n as f64).sqrt() as usize).clamp(5, 100).min(n);
    let width = max_s / bins as f64;

    let mut cells: Vec<(f64, f64)> = Vec::with_capacity(bins);
    let mut observed = vec![0.0_f64; bins];
    for &s in spacings {
        let k = (s / width).floor() as usize;
        observed[k.min(bins - 1)] += 1.0;
    }
    for (i, &o) in observed.iter().enumerate() {
        let lo = i as f64 * width;
        let hi = (i + 1) as f64 * width;
        let prob = (-lo).exp() - if i + 1 == bins { 0.0 } else { (-hi).exp() };
        cells.push((o, prob * n as f64));
    }

    // Merge leftward until every non-terminal cell expects at least 5 counts.
    let mut merged: Vec<(f64, f64)> = Vec::with_capacity(cells.len());
    for cell in cells {
        merged.push(cell);
        while merged.len() >= 2 && merged[merged.len() - 2].1 < 5.0 {
            let last = merged.pop().expect("len >= 2");
            let prev = merged.last_mut().expect("len >= 1");
            prev.0 += last.0;
            prev.1 += last.1;
        }
    }

    let chi2: f64 = merged
        .iter()
        .filter(|(_, e)| *e > 0.0)
        .map(|(o, e)| (o - e).powi(2) / e)
        .sum();
    let df = merged.len().saturating_sub(2);
    if df == 0 {
        return Err("not enough bins after merging".into());
    }
    let p = gammq(df as f64 / 2.0, chi2 / 2.0);
    Ok((chi2, df, p))
}

/// Full spectral test on one candidate adjacency.
pub fn spectral_test(adj: &DMatrix<f64>, threshold: f64) -> Result<Probe, MenapError> {
    let eig = adj.clone().symmetric_eigen();
    let mut evs: Vec<f64> = eig
        .eigenvalues
        .iter()
        .copied()
        .filter(|x| x.abs() > ZERO_EIG)
        .collect();
    evs.sort_by(f64::total_cmp);

    if evs.len() < MIN_LEVELS {
        return Err(MenapError::TooFewEigenvalues {
            threshold,
            levels: evs.len(),
        });
    }
    if !evs.iter().all(|x| x.is_finite()) {
        return Err(MenapError::Degenerate {
            threshold,
            what: "non-finite eigenvalue".into(),
        });
    }

    let unfolded = unfold(&evs);
    let d = spacings(&unfolded);
    if d.len() < MIN_LEVELS {
        return Err(MenapError::Degenerate {
            threshold,
            what: "too many repeated unfolded levels".into(),
        });
    }
    let (chi2, df, p) =
        chi_square_poisson(&d).map_err(|what| MenapError::Degenerate { threshold, what })?;

    Ok(Probe {
        threshold,
        nodes: adj.nrows(),
        edges: 0,
        chi2,
        df,
        p,
        accepted: false,
    })
}

#[derive(Debug, Clone)]
pub struct ScanConfig {
    pub min: f64,
    pub max: f64,
    pub coarse: f64,
    pub fine: f64,
    pub alpha: f64,
}

impl Default for ScanConfig {
    fn default() -> Self {
        Self {
            min: SEED,
            max: SCAN_MAX,
            coarse: COARSE_STEP,
            fine: FINE_STEP,
            alpha: ALPHA,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Outcome {
    pub threshold: f64,
    pub network: Network,
    pub probe: Probe,
    /// True when the scan accepted at its lower bound, so the transition may
    /// lie below the scanned range.
    pub at_lower_bound: bool,
    /// True when no coarse threshold was compatible with Poisson and the whole
    /// range had to be swept at the fine step (the accept window is narrower
    /// than the coarse step).
    pub full_sweep: bool,
    /// Every evaluated threshold, in order.
    pub curve: Vec<Probe>,
}

fn evaluate(sim: &DMatrix<f64>, t: f64, alpha: f64) -> Result<(Network, Probe), MenapError> {
    let net = build_network(sim, t);
    if net.is_empty() {
        return Err(MenapError::EmptyNetwork { threshold: t });
    }
    let mut probe = spectral_test(&net.adj, t)?;
    probe.nodes = net.len();
    probe.edges = net.edges();
    probe.accepted = probe.p > alpha;
    Ok((net, probe))
}

/// Two-stage scan, as in the paper: coarse steps until the test stops
/// rejecting, then fine steps inside that bracket. The first accepted
/// threshold wins.
///
/// The coarse pass assumes the statistic moves monotonically with the
/// threshold. It often does not (the spectrum of a sparse random graph can
/// look Poisson-like at one threshold and Wigner-like again at a denser one),
/// so when no coarse candidate is accepted the scan falls back to a full sweep
/// at the fine step. That fallback is what catches accept windows narrower
/// than the coarse grid.
///
/// Empty and degenerate candidates are terminal, not skipped: raising the
/// threshold removes edges and levels monotonically, so a later candidate
/// cannot recover.
pub fn scan(sim: &DMatrix<f64>, cfg: &ScanConfig) -> Result<Outcome, MenapError> {
    if cfg.fine <= 0.0 || cfg.coarse <= 0.0 || cfg.max <= cfg.min {
        return Err(MenapError::Degenerate {
            threshold: cfg.min,
            what: "invalid scan bounds or steps".into(),
        });
    }

    let mut curve: Vec<Probe> = Vec::new();
    let mut seen: std::collections::HashSet<i64> = std::collections::HashSet::new();
    let mut last_p = f64::NAN;
    let mut failure: Option<String> = None;

    // --- pass 1: coarse sweep ------------------------------------------------
    let coarse_steps = ((cfg.max - cfg.min) / cfg.coarse).floor() as usize;
    let mut previous = cfg.min;
    let mut bracket: Option<(f64, f64, Network, Probe)> = None;

    for k in 0..=coarse_steps {
        let t = round6(cfg.min + k as f64 * cfg.coarse);
        seen.insert(key(t));
        match evaluate(sim, t, cfg.alpha) {
            Ok((net, probe)) => {
                last_p = probe.p;
                let accepted = probe.accepted;
                curve.push(probe.clone());
                if accepted {
                    if k == 0 {
                        return Ok(Outcome {
                            threshold: t,
                            network: net,
                            probe,
                            at_lower_bound: true,
                            full_sweep: false,
                            curve: sorted(curve),
                        });
                    }
                    bracket = Some((previous, t, net, probe));
                    break;
                }
                previous = t;
            }
            Err(e) => {
                failure = Some(e.to_string());
                break;
            }
        }
    }

    // --- pass 2a: refine inside the coarse bracket --------------------------
    if let Some((lo, hi, net, probe)) = bracket {
        let steps = ((hi - lo) / cfg.fine).round() as usize;
        for k in 1..steps {
            let t = round6(lo + k as f64 * cfg.fine);
            if !seen.insert(key(t)) {
                continue;
            }
            match evaluate(sim, t, cfg.alpha) {
                Ok((net, probe)) => {
                    let accepted = probe.accepted;
                    curve.push(probe.clone());
                    if accepted {
                        return Ok(Outcome {
                            threshold: t,
                            network: net,
                            probe,
                            at_lower_bound: false,
                            full_sweep: false,
                            curve: sorted(curve),
                        });
                    }
                }
                // A candidate inside the bracket cannot be sparser than the
                // bracket's upper end, which passed; stop refining and return
                // that end.
                Err(_) => break,
            }
        }
        // The bracket's upper end is known to pass and nothing below it in the
        // bracket did, so it is the refined answer.
        return Ok(Outcome {
            threshold: hi,
            network: net,
            probe,
            at_lower_bound: false,
            full_sweep: false,
            curve: sorted(curve),
        });
    }

    // --- pass 2b: full fine sweep -------------------------------------------
    if failure.is_none() {
        let steps = ((cfg.max - cfg.min) / cfg.fine).floor() as usize;
        for k in 1..=steps {
            let t = round6(cfg.min + k as f64 * cfg.fine);
            if !seen.insert(key(t)) {
                continue;
            }
            match evaluate(sim, t, cfg.alpha) {
                Ok((net, probe)) => {
                    last_p = probe.p;
                    let accepted = probe.accepted;
                    curve.push(probe.clone());
                    if accepted {
                        return Ok(Outcome {
                            threshold: t,
                            network: net,
                            probe,
                            at_lower_bound: false,
                            full_sweep: true,
                            curve: sorted(curve),
                        });
                    }
                }
                Err(e) => {
                    failure = Some(e.to_string());
                    break;
                }
            }
        }
    }

    Err(MenapError::NoThreshold {
        max: cfg.max,
        last_p,
        reason: failure.unwrap_or_else(|| "the test rejected every candidate".into()),
    })
}

/// Quantise a threshold to 1e-6 so grid points compare exactly and sort in
/// order (0.30 + 56 * 0.01 must print and compare as 0.86).
fn round6(x: f64) -> f64 {
    (x * 1e6).round() / 1e6
}

fn key(t: f64) -> i64 {
    (t * 1e6).round() as i64
}

/// Curves are assembled from two passes; present them in threshold order.
fn sorted(mut curve: Vec<Probe>) -> Vec<Probe> {
    curve.sort_by(|a, b| a.threshold.total_cmp(&b.threshold));
    curve
}

// ---------------------------------------------------------------------------
// numerics
// ---------------------------------------------------------------------------

/// Complementary error function: Numerical Recipes' rational approximation,
/// fractional error below 1.2e-7 across the whole range.
pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07
                                + t * (-1.135_203_98
                                    + t * (1.488_515_87
                                        + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}

/// Error function, `1 - erfc(x)`.
pub fn erf(x: f64) -> f64 {
    1.0 - erfc(x)
}

/// Standard normal CDF, `Phi(x) = erfc(-x / sqrt(2)) / 2`.
pub fn normal_cdf(x: f64) -> f64 {
    0.5 * erfc(-x / std::f64::consts::SQRT_2)
}

/// Regularised upper incomplete gamma `Q(a, x)`, i.e. the chi-square tail
/// probability when `a = df/2` and `x = chi2/2`.
pub fn gammq(a: f64, x: f64) -> f64 {
    if a.is_nan() || a <= 0.0 || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return 1.0;
    }
    if x < a + 1.0 {
        1.0 - gser(a, x)
    } else {
        gcf(a, x)
    }
}

fn gser(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..1000 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-15 {
            break;
        }
    }
    sum * (-x + a * x.ln() - gln).exp()
}

fn gcf(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    const TINY: f64 = 1e-300;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / TINY;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..1000 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < TINY {
            d = TINY;
        }
        c = b + an / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-15 {
            break;
        }
    }
    (-x + a * x.ln() - gln).exp() * h
}

fn ln_gamma(x: f64) -> f64 {
    const COF: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection; only used defensively, our a is always >= 0.5.
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = COF[0];
    let t = x + 7.5;
    for (i, &c) in COF.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_file(name: &str, body: &str) -> std::path::PathBuf {
        let dir =
            std::env::temp_dir().join(format!("rusty-menap-test-{}-{}", std::process::id(), name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("table.txt");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        path
    }

    #[test]
    fn loader_accepts_crlf_and_missing_cells() {
        let p = temp_file("crlf", "ID\tA\tB\tC\r\nOTU1\t1\t\t3\r\nOTU2\t4\t5\t6\r\n");
        let ds = load_otu_table(&p).unwrap();
        assert_eq!(ds.ids, vec!["OTU1", "OTU2"]);
        assert_eq!(ds.values.shape(), (2, 3));
        assert_eq!(ds.values[(0, 1)], 0.0);
        assert_eq!(ds.values[(1, 2)], 6.0);
    }

    #[test]
    fn loader_skips_blank_lines() {
        let p = temp_file("blank", "ID\tA\tB\n\nOTU1\t1\t2\n\nOTU2\t3\t4\n");
        let ds = load_otu_table(&p).unwrap();
        assert_eq!(ds.ids.len(), 2);
        assert_eq!(ds.values.shape(), (2, 2));
    }

    #[test]
    fn loader_rejects_ragged_rows() {
        let p = temp_file("ragged", "ID\tA\tB\nOTU1\t1\t2\nOTU2\t3\t4\t5\n");
        match load_otu_table(&p) {
            Err(MenapError::Ragged {
                line,
                expected,
                found,
            }) => {
                assert_eq!((line, expected, found), (3, 3, 4));
            }
            other => panic!("expected Ragged, got {other:?}"),
        }
    }

    #[test]
    fn loader_rejects_non_numeric_cells() {
        let p = temp_file("bad", "ID\tA\tB\nOTU1\t1\tx\n");
        assert!(matches!(
            load_otu_table(&p),
            Err(MenapError::BadNumber { column: 3, .. })
        ));
    }

    #[test]
    fn loader_rejects_header_only() {
        let p = temp_file("header", "ID\tA\tB\n");
        assert!(matches!(load_otu_table(&p), Err(MenapError::NoData)));
    }

    #[test]
    fn standardize_produces_zero_mean_unit_sd() {
        let mut v = DMatrix::zeros(2, 4);
        for (i, x) in [1.0, 2.0, 3.0, 4.0].iter().enumerate() {
            v[(0, i)] = *x;
            v[(1, i)] = *x * 2.0 + 1.0;
        }
        let (z, constant) = standardize(&v).unwrap();
        assert!(constant.is_empty());
        for r in 0..2 {
            let mean: f64 = z.row(r).iter().sum::<f64>() / 4.0;
            let sd = (z.row(r).iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 3.0).sqrt();
            assert!(mean.abs() < 1e-12);
            assert!((sd - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn standardize_flags_constant_rows() {
        let mut v = DMatrix::zeros(2, 3);
        for c in 0..3 {
            v[(0, c)] = 7.0;
            v[(1, c)] = c as f64;
        }
        let (_, constant) = standardize(&v).unwrap();
        assert_eq!(constant, vec![0]);
    }

    #[test]
    fn similarity_is_symmetric_bounded_and_unit_diagonal() {
        let mut v = DMatrix::zeros(4, 6);
        for r in 0..4 {
            for c in 0..6 {
                v[(r, c)] = ((r * 7 + c * 3) % 5) as f64 + r as f64 * 0.1;
            }
        }
        let (z, _) = standardize(&v).unwrap();
        let s = similarity(&z);
        for i in 0..4 {
            assert!((s[(i, i)] - 1.0).abs() < 1e-15);
            for j in 0..4 {
                assert!((s[(i, j)] - s[(j, i)]).abs() < 1e-12);
                assert!((0.0..=1.0).contains(&s[(i, j)]));
            }
        }
    }

    #[test]
    fn similarity_makes_anticorrelation_positive() {
        let mut v = DMatrix::zeros(2, 4);
        for (c, x) in [1.0, 2.0, 3.0, 4.0].iter().enumerate() {
            v[(0, c)] = *x;
            v[(1, c)] = -*x;
        }
        let (z, _) = standardize(&v).unwrap();
        let s = similarity(&z);
        assert!((s[(0, 1)] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn network_drops_isolates_and_zeroes_the_diagonal() {
        let mut sim = DMatrix::from_element(4, 4, 0.9);
        // node 3 is below threshold with every other node
        for j in 0..4 {
            sim[(3, j)] = 0.1;
            sim[(j, 3)] = 0.1;
        }
        for i in 0..4 {
            sim[(i, i)] = 1.0;
        }
        let net = build_network(&sim, 0.5);
        assert_eq!(net.nodes, vec![0, 1, 2]);
        assert_eq!(net.len(), 3);
        for i in 0..3 {
            assert_eq!(net.adj[(i, i)], 0.0);
        }
        assert_eq!(net.edges(), 3);
        assert!((net.mean_degree() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn unfolding_matches_the_naive_full_sum() {
        // Deterministic pseudo-spectrum with a realistic gap structure.
        let mut evs: Vec<f64> = (0..60)
            .map(|i| (i as f64 * 0.37).sin() * 2.0 + i as f64 * 0.1)
            .collect();
        evs.sort_by(f64::total_cmp);
        let fast = unfold(&evs);
        // Naive O(n^2) reference: same kernel, no window and no (k+1)
        // decomposition, so this checks exactly the algebra the fast path
        // relies on.
        let n = evs.len();
        let sigma = SIGMA_LEVELS * (evs[n - 1] - evs[0]) / n as f64;
        for (k, &e) in evs.iter().enumerate() {
            let naive: f64 = evs.iter().map(|&ej| normal_cdf((e - ej) / sigma)).sum();
            assert!(
                // The kernel itself is accurate to ~1.2e-7 (NR erfc), which is
                // the floor for this comparison.
                (fast[k] - naive).abs() < 1e-7,
                "level {k}: {} vs {naive}",
                fast[k]
            );
        }
    }

    #[test]
    fn erfc_is_accurate_at_the_textbook_points() {
        assert!((erf(0.0)).abs() < 1e-6);
        assert!((erf(1.0) - 0.842_700_792_948_629_5).abs() < 2e-7);
        assert!((erfc(1.0) - 0.157_299_207_050_285_13).abs() < 2e-8);
        assert!((normal_cdf(0.0) - 0.5).abs() < 1e-7);
        assert!((normal_cdf(1.959_963_984_540_054) - 0.975).abs() < 1e-7);
    }

    #[test]
    fn spacings_have_unit_mean() {
        let levels: Vec<f64> = (0..50).map(|i| i as f64 * i as f64 * 0.01).collect();
        let d = spacings(&levels);
        let mean = d.iter().sum::<f64>() / d.len() as f64;
        assert!((mean - 1.0).abs() < 1e-12);
    }

    /// Small deterministic LCG so the statistical tests do not need a RNG
    /// dependency (and cannot flake on a different runtime).
    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64)
    }

    #[test]
    fn chi_square_accepts_poisson_spacings() {
        let mut seed = 42;
        let d: Vec<f64> = (0..900).map(|_| -(1.0 - lcg(&mut seed)).ln()).collect();
        let (_, df, p) = chi_square_poisson(&d).unwrap();
        assert!(df > 0);
        assert!(p > 0.01, "Poisson sample rejected at p = {p}");
    }

    #[test]
    fn chi_square_rejects_wigner_dyson_spacings() {
        // GOE Wigner surmise: F(s) = 1 - exp(-pi s^2 / 4).
        let mut seed = 7;
        let d: Vec<f64> = (0..900)
            .map(|_| (-4.0 * (1.0 - lcg(&mut seed)).ln() / std::f64::consts::PI).sqrt())
            .collect();
        let (_, _, p) = chi_square_poisson(&d).unwrap();
        assert!(p < 0.01, "Wigner-Dyson sample accepted at p = {p}");
    }

    #[test]
    fn gammq_matches_closed_forms() {
        assert!((gammq(1.0, 2.0) - (-2.0_f64).exp()).abs() < 1e-12);
        // Q(1/2, x) = erfc(sqrt(x))
        assert!((gammq(0.5, 1.0) - 0.157_299_207_050_285_13).abs() < 2e-8);
        // chi-square df=30 at the mean has p = 0.4657
        let p = gammq(15.0, 15.0);
        assert!((p - 0.4657).abs() < 1e-3, "p = {p}");
    }

    #[test]
    fn scan_returns_a_threshold_and_a_consistent_network() {
        // 60 OTUs x 24 samples in 3 planted blocks plus deterministic noise.
        let (n, m) = (60, 24);
        let mut v = DMatrix::zeros(n, m);
        let mut seed = 12345;
        for i in 0..n {
            let block = i / 20;
            for c in 0..m {
                let base = (c as f64 * 0.7 + block as f64).sin();
                v[(i, c)] = base + 0.35 * lcg(&mut seed);
            }
        }
        let (z, _) = standardize(&v).unwrap();
        let sim = similarity(&z);
        let out = scan(&sim, &ScanConfig::default());
        if let Ok(out) = out {
            assert!(out.threshold >= SEED && out.threshold < SCAN_MAX);
            assert!(out.network.len() <= n);
            assert_eq!(out.network.edges(), out.probe.edges);
            assert_eq!(out.network.len(), out.probe.nodes);
            for i in 0..out.network.len() {
                assert_eq!(out.network.adj[(i, i)], 0.0);
            }
        }
    }
}
