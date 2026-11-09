# Rusty-MENAP

Rust implementation of the random-matrix-theory (RMT) threshold search used by
the Molecular Ecological Network Analysis pipeline (MENAP) of
Deng et al. (2012), *Molecular ecological network analyses*, BMC Bioinformatics
13:113 ([doi:10.1186/1471-2105-13-113](https://doi.org/10.1186/1471-2105-13-113)).

Given an OTU abundance table (rows = OTUs, columns = samples), the tool
standardises each OTU, computes the absolute Pearson correlation matrix, and
scans candidate thresholds until the spectrum of the thresholded adjacency
matrix is statistically compatible with a random (Poisson) one. It writes the
resulting edge list, the fitted statistics, and the full scan curve.

```
input table ──► standardise ──► S = |R| ──► threshold t ──► A(t) ──► eigenvalues
                                                  │                     │
                                                  │              unfold + NNSD
                                                  │                     │
                                                  └──── scan ◄──── χ² vs Poisson
                                                          │
                                             first t that does not reject
```

## Build and run

```bash
cargo build --release
./target/release/rusty-menap                      # data/Control.txt, out/
./target/release/rusty-menap -i data/warming.txt -o out-warming
cargo run --release -- --st 0.52                  # fixed threshold, no scan
```

Requires Rust 1.88 or newer (edition 2024). Pixi users: `pixi run build`,
`pixi run start`.

## Options

| flag | default | meaning |
|---|---|---|
| `-i, --input` | `data/Control.txt` | OTU table, tab-separated, header row required |
| `-o, --out` | `out` | output directory |
| `--alpha` | `0.01` | significance level of the χ² goodness-of-fit test |
| `--st` | – | use a fixed threshold instead of scanning |
| `--min`, `--max` | `0.30`, `0.99` | scan bounds |
| `--coarse`, `--fine` | `0.10`, `0.01` | scan steps |
| `--no-curve` | off | do not write `scan_curve.tsv` |
| `--no-write` | off | analyse without writing anything |
| `-q, --quiet` | off | print only the final summary |

Exit codes: `0` success, `1` analysis failure (no threshold, degenerate
spectrum, malformed input), `2` usage error.

## Outputs

- `network.tsv` — edge list `source<TAB>target<TAB>weight`, one row per
  undirected edge, weights are the correlation magnitudes that survived.
- `summary.json` — input, mode, α, chosen threshold, node/edge counts, mean
  degree, and the χ², df and p of the accepted candidate.
- `scan_curve.tsv` — one row per evaluated threshold: nodes, edges, χ², df, p,
  decision. Useful for seeing whether the transition is clean or marginal.

## Method, and where it is a choice rather than the paper

Steps (a)–(h) of the paper's Methods are implemented as written: row-wise
standardisation, `S = |R|`, weighted adjacency `a_ij = s_ij` above the
threshold (Eq. 3), isolated nodes dropped so that `p ≤ n`, eigenvalues of the
adjacency, unfolding, nearest-neighbour spacing distribution (NNSD), χ²
goodness-of-fit against `exp(-d)`.

The paper leaves several details open. This implementation's choices are:

- **Unfolding** by Gaussian broadening with `σ = 5 × mean level spacing`
  (paper step f allows a "local density average"; this is one).
- **Binning**: `round(√n)` equal-width bins up to the largest observed spacing,
  the last bin open-ended, adjacent bins merged until each expectation reaches
  5 — so no single cell dominates the statistic.
- **Degrees of freedom** `bins - 2`, because the mean spacing is estimated
  from the same data. At α = 0.01 the acceptance region is the upper tail,
  `p > 0.01`.
- **Exact null space removed**: eigenvalues within `1e-9` of zero (one per
  disconnected component) are excluded from the bulk.
- **Scan strategy**: coarse sweep, then a fine sweep inside the first coarse
  accept bracket. If no coarse candidate passes — the statistic is *not*
  monotonic in the threshold, and accept windows narrower than 0.10 do occur —
  the whole range is re-swept at the fine step (several times slower, and
  reported as such).
- **Threshold resolution** is the fine step; the reported value is the first
  grid point that passes, not an interpolated crossing.

## What the test does and does not prove

A pass means "this adjacency spectrum is not distinguishable from a Poisson
one". It does **not** mean the network is biologically meaningful: on a
correlation-free null of the same shape (every OTU's values permuted across
samples) the same scan still returns a threshold, because sparsity alone drives
the spectrum towards Poisson. Treat the threshold as a scan statistic that
needs a null model and a structural check (modularity against degree-matched
random graphs, module stability) before it is interpreted as signal. The
reference pipeline itself abandons RMT for a fixed threshold when the test does
not accept on its data.

Runtime on the shipped control table (987 OTUs × 14 samples): about 8 s when
the coarse scan brackets the transition, about 35 s when the full fine sweep is
needed. Memory is O(n²) f64 (≈8 MB at n = 1000, ≈800 MB at n = 10 000).

## Example data

`data/Control.txt` and `data/warming.txt` are the control and warming OTU
tables distributed with the paper. `ressources/MENA.pdf` is the paper itself.

## Layout

```
src/lib.rs      loader, standardisation, similarity, adjacency, RMT test, scan
src/main.rs     CLI, output writers
tests/cli.rs    end-to-end behaviour of the binary
examples/probe.rs  print the spectral probe at fixed thresholds
```

## License

GPL-3.0-or-later; see `LICENSE`.
