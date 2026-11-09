# Rusty-MENAP — review and fix status

**Date:** 2026-11-09 · **Scope:** the 399-line prototype at `src/main.rs` (commit `257d04a`,
whole implementation uncommitted) and its rewrite in `0.2.0`.

Reference: Deng et al. 2012, *Molecular ecological network analyses*, BMC Bioinformatics
13:113 (`ressources/MENA.pdf`), Methods steps (a)–(h).

Every number below was produced by running the code on the shipped
`data/Control.txt` (987 OTUs × 14 samples), not estimated.

## Baseline (0.1.0)

```
$ ./target/release/Rusty-MENAP
Optimal Threshold found: 0.84
Final Network: 987 Nodes, 2842 Edges
```

## Findings and status

| # | Finding | Severity | Status |
|---|---|---|---|
| 1 | Reported edge count inflated by n/2: the 1.0 diagonal was counted as a positive entry and halved (2842 reported, 2349 true) | high | fixed |
| 2 | Isolated nodes never pruned, contradicting step (d); 987 nodes reported while only 786 had an edge | high | fixed |
| 3 | 1.0 diagonal made every row non-zero, so the "network disintegrated" guard was unreachable | high | fixed |
| 4 | Paper-faithful adjacency (zero diagonal, pruning) moves the transition 0.84 → 0.86, a 27 % change in edge count | high | fixed |
| 5 | The criterion does not separate signal from noise: a correlation-free null (each OTU permuted across samples) returned 0.83 against 0.84 on the real table | high | documented, not fixable in code |
| 6 | Acceptance was one-sided through a normal approximation of χ² and accepted "too good" fits; no exact p-value, df = bins − 1 with the mean spacing estimated | medium | fixed |
| 7 | The statistic is non-monotonic in the threshold, so the coarse-then-fine scan could skip a narrow accept window | medium | fixed (full fine-sweep fallback) |
| 8 | Threshold sensitive to arbitrary σ and bin count (σ 1–8 → 0.84–0.86; bins 10–100 → 0.84–0.90) | medium | documented; σ and binning are explicit choices |
| 9 | "Too few eigenvalues" raised the threshold, making the graph sparser and the problem worse | medium | fixed (terminal error) |
| 10 | Refine step was a no-op when the coarse scan passed at the seed: it decremented, clamped to 0.30 and re-tested the same point | medium | fixed (two-stage scan, `at_lower_bound` flag) |
| 11 | Loader fabricated an OTU for a blank line, panicked on a row longer than the header, and turned unparsable cells into 0.0 | medium | fixed (strict, line-numbered errors) |
| 12 | `cargo clippy` exited 101 (`0./0.`), `cargo fmt --check` failed | medium | fixed (0 warnings) |
| 13 | Unused `rand` dependency and three unused imports; duplicated section comment | low | fixed |
| 14 | No README, no tests (`cargo test`: 0), no output files, hardcoded input path, `warming.txt` unused | low | fixed |
| 15 | `max = 0.99` was never probed (`while threshold < 0.99`) | low | fixed |

## What changed

- **`src/lib.rs`** (new): loader, standardisation, similarity, adjacency,
  RMT test and scan as a testable library. Errors carry line numbers.
- **`src/main.rs`**: CLI (`--input/--out/--alpha/--st/--min/--max/--coarse/--fine`),
  edge-list, JSON summary and scan-curve writers.
- **`build_network`**: masks the similarity matrix strictly above the diagonal,
  drops OTUs with no surviving edge, and zeroes the diagonal.
- **`Network::edges`**: counts only the upper triangle.
- **`chi_square_poisson`**: exact χ² tail via the regularised incomplete gamma
  (series + continued fraction), bins merged to expectation ≥ 5, df = bins − 2.
- **`unfold`**: same Gaussian-broadened cumulative density, but accumulated with
  `erfc` on both tails (no `1 − erf` cancellation) and truncated at 8σ; equals
  the naive O(n²) sum to the kernel's 1.2e-7.
- **`scan`**: coarse sweep → fine sweep in the bracket → full fine sweep when no
  coarse candidate passes; errors are terminal instead of raising the threshold.
- **`tests/cli.rs`**: end-to-end checks of outputs and exit codes.
- **`examples/probe.rs`**: prints the spectral probe at chosen thresholds.

## Verification

```
$ cargo test
17 unit tests + 4 CLI tests, all passing

$ cargo clippy --all-targets --release
0 warnings            $ cargo fmt --check → clean

$ ./target/release/rusty-menap
threshold 0.86  (rmt mode)
network   695 nodes / 1713 edges, mean degree 4.93, 292 isolates dropped
fit       chi2 = 22.76, df = 15, p = 0.0894 > alpha = 0.01
```

The edge count now equals the number of lines in `out/network.tsv` (1713), the
node count equals the connected OTUs at that threshold, and the reported
threshold matches an independent NumPy reimplementation of the corrected
procedure (`0.86`, 695 nodes, 1713 edges).

Reproduction of finding 5 on the corrected build: permuting every OTU's values
across samples yields `threshold 0.83`, `p = 0.0113` — still accepted at
α = 0.01. At α = 0.05 the null is rejected (0.0113) while the real table passes
(0.0894), so the default α = 0.01 follows the paper but does not by itself
separate this control table from its null.

## Not fixed

- **No null model.** The scan statistic alone cannot establish that a threshold
  is meaningful. Structural checks (modularity against degree-matched random
  graphs, bootstrap module stability) are outside this tool's scope.
- **Lower-tail rejection.** A χ² far below its expectation (spacings more
  regular than Poisson) is still accepted; the paper defines only the upper
  tail, so adding a two-sided rule would invent method.
- **Scaling.** Dense O(n³) eigendecomposition and O(n²) f64 storage: about 8 s
  per probe at n = 1000 and 800 MB per matrix at n = 10 000. No sparse path.
- **Cross-validation with the MENAP R implementation** (not available here):
  numeric agreement is unverified.
