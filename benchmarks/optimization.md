# glyrex performance optimization — 2026-10-24

Input structure parsing is unchanged: character inputs still call the same `glyparse::auto_parse()` with the same arguments and inputs. No dependency was upgraded. Improvements are in pattern compilation, matching and output construction.

## Paired measurement

Old and new packages were installed in separate temporary libraries and measured sequentially in fresh R processes with identical dependencies. Seven timed batches per case; median wall time includes garbage collection during measurement. Package loading, input preparation and precompilation are excluded where indicated. There is no cross-call cache of query results; only immutable residue metadata is cached across calls.

The shared pattern is `(?<terminal>Hex)-HexNAc`. The 100-input workload repeats 26 fixture glycans (2–15 residues); the unique control uses each of those glycans once. These fixtures are not a population-weighted biological dataset. Structure extraction includes construction of `glyrepr_structure` outputs.

## 100 structures, precompiled pattern

| Case | Before (ms) | After (ms) | Speedup |
|---|---:|---:|---:|
| rex_count | 66.000 | 1.719 | 38.40x |
| rex_detect | 69.000 | 1.719 | 40.15x |
| rex_extract | 92.000 | 3.125 | 29.44x |
| rex_extract_all | 140.000 | 5.375 | 26.05x |
| rex_locate | 71.000 | 2.063 | 34.42x |
| rex_locate_all | 72.000 | 3.375 | 21.33x |
| rex_match | 120.000 | 3.687 | 32.54x |
| rex_match_all | 177.000 | 4.937 | 35.85x |
| rex_subset | 64.000 | 1.844 | 34.71x |
| rex_which | 65.000 | 1.750 | 37.14x |

## 100 strings, string pattern (including unchanged parsing)

| Case | Before (ms) | After (ms) | Speedup |
|---|---:|---:|---:|
| rex_count | 95.000 | 22.000 | 4.32x |
| rex_detect | 95.000 | 22.000 | 4.32x |
| rex_extract | 112.000 | 23.000 | 4.87x |
| rex_extract_all | 142.000 | 23.500 | 6.04x |
| rex_locate | 96.000 | 22.000 | 4.36x |
| rex_locate_all | 93.000 | 25.000 | 3.72x |
| rex_match | 146.000 | 23.750 | 6.15x |
| rex_match_all | 194.000 | 25.500 | 7.61x |
| rex_subset | 86.000 | 22.500 | 3.82x |
| rex_which | 88.000 | 22.250 | 3.96x |

## 26 distinct structures, precompiled pattern

| Case | Before (ms) | After (ms) | Speedup |
|---|---:|---:|---:|
| rex_count | 16.750 | 1.719 | 9.75x |
| rex_detect | 16.750 | 1.719 | 9.75x |
| rex_extract | 24.000 | 3.062 | 7.84x |
| rex_extract_all | 35.500 | 5.313 | 6.68x |
| rex_locate | 18.750 | 1.875 | 10.00x |
| rex_locate_all | 18.250 | 3.375 | 5.41x |
| rex_match | 31.500 | 3.563 | 8.84x |
| rex_match_all | 44.000 | 4.812 | 9.14x |
| rex_subset | 16.750 | 1.750 | 9.57x |
| rex_which | 16.750 | 1.656 | 10.11x |

## Pattern compilation

| Case | Before (ms) | After (ms) | Speedup |
|---|---:|---:|---:|
| absence | 1.719 | 0.168 | 10.23x |
| assertion | 3.437 | 0.371 | 9.26x |
| branch | 1.812 | 0.266 | 6.82x |
| capture | 1.203 | 0.152 | 7.90x |
| linkage | 1.172 | 0.148 | 7.89x |
| no_match | 0.578 | 0.063 | 9.11x |
| repeat_path | 0.609 | 0.102 | 6.00x |
| simple | 0.578 | 0.062 | 9.25x |

## Detection across pattern types (100 structures)

| Case | Before (ms) | After (ms) | Speedup |
|---|---:|---:|---:|
| absence | 96.000 | 6.000 | 16.00x |
| assertion | 99.000 | 4.812 | 20.57x |
| branch | 65.000 | 4.250 | 15.29x |
| capture | 64.000 | 1.719 | 37.24x |
| linkage | 55.000 | 2.844 | 19.34x |
| no_match | 48.500 | 1.141 | 42.52x |
| repeat_path | 60.000 | 2.719 | 22.07x |
| simple | 57.000 | 1.719 | 33.16x |

## Detection across batch sizes

| Case | Before (ms) | After (ms) | Speedup |
|---|---:|---:|---:|
| 1 structure_compiled | 0.727 | 0.250 | 2.91x |
| 1 text | 2.937 | 1.406 | 2.09x |
| 10 structure_compiled | 7.000 | 0.797 | 8.78x |
| 10 text | 16.750 | 9.125 | 1.84x |
| 100 structure_compiled | 65.000 | 1.812 | 35.86x |
| 100 text | 86.000 | 22.500 | 3.82x |
| 1000 structure_compiled | 649.000 | 2.500 | 259.60x |
| 1000 text | 670.000 | 23.000 | 29.13x |

## Implementation and verification

- Each unique structure-pattern pair is matched once per call, then results are restored to the original order and names.
- Fixed paths use vectorized matching; general backtracking retains its semantics and reuses atom masks within a query.
- Raw igraph inputs retain full validation. Existing structure vectors reuse their validated graphs; floating metadata remains rejected.
- Fragments are serialized in batches; ambiguous branch-order ties use the original graph extraction path to preserve exact strings.
- 790 baseline-versus-optimized result objects agree exactly after normalizing structure objects to their typed, named IUPAC values. They cover all 10 query functions, 39 patterns, 26 fixture glycans, 30 seeded random glycans (2–8 residues), missing values, and permuted raw graphs.
- Package tests: 485 passed. R CMD check: 0 errors, 0 warnings, 0 notes (`--no-manual`, `_R_CHECK_SYSTEM_CLOCK_=FALSE` to disable unavailable remote clock verification).
- Fixed-path tests include 192 assertions comparing the fast path with the backtracking engine. A discovered ambiguous-branch serialization difference has a dedicated regression test.
- All 39 compiled syntax trees are identical to the baseline after the final token extraction change.

## Reproduce

Install each source revision into a separate R library, then run:

```sh
Rscript benchmarks/functions.R benchmarks/results-before-optimization /path/to/before-library BEFORE_REVISION
Rscript benchmarks/functions.R benchmarks/results-after-optimization /path/to/after-library AFTER_REVISION
Rscript benchmarks/compare.R
Rscript benchmarks/parity.R /path/to/before-library /tmp/before.rds
Rscript benchmarks/parity.R /path/to/after-library /tmp/after.rds
Rscript -e 'stopifnot(identical(readRDS("/tmp/before.rds"), readRDS("/tmp/after.rds")))'
```

Raw batch timings, corpus copies, dependency versions and source revisions are retained in `results-before-optimization/` and `results-after-optimization/`. `comparison.csv` contains all 64 paired cases. The original first-run baseline remains in `results-2026-10-24/`.
