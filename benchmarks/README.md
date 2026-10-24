# glyrex function benchmark — 2026-10-24

This is the historical baseline. See [the optimization comparison](optimization.md)
for paired before/after results. New runs default to `results-current/` so that
the original baseline is preserved.

Run from the package root with `Rscript benchmarks/functions.R`.
The script loads the working-tree package using `devtools::load_all()`.

## Method

- All 11 exported functions; 54 benchmark cases with seven timed batches each.
- Input corpus: unique glycan strings from `tests/testthat/fixtures/glycowork.R`, saved in `results-2026-10-24/corpus.csv`.
- API comparison: repeat this corpus to 100 glycans, using `(?<terminal>Hex)-HexNAc` for every query function.
- Additional cases: eight pattern types and detection at 1, 10, 100, and 1,000 glycans.
- Warm-session wall time, including garbage collection during the timed batch. Each case is warmed and calibrated to at least 50 ms per batch; explicit GC occurs before each batch outside timing. Report the median per-call time across seven batches.
- Input creation and precompilation are outside timing where applicable. String-input calls include parsing and, for text patterns, compilation. Package loading and cold-start cost are excluded.
- Repeated fixtures are a reproducible workload, not a representative biological dataset. Parsing can benefit from duplicate inputs. Scaling cases below the corpus size use its first entries, so their composition differs.
- Structure input changes `extract`/`extract_all` return types to structure vectors, so their comparison includes that output conversion.
- These sequential desktop timings are approximate, not randomized controlled measurements. The short test-suite run overlapped the initial API cases; interpret small differences as noise.

## Results

Median milliseconds per call on 100 glycans:

| Function | String + string pattern | String + compiled pattern | Structure + compiled pattern |
|---|---:|---:|---:|
| rex_detect | 97 | 103 | 73 |
| rex_count | 92 | 93 | 73 |
| rex_which | 92 | 89 | 68 |
| rex_subset | 100 | 90 | 69 |
| rex_locate | 100 | 99 | 79 |
| rex_locate_all | 98 | 96 | 76 |
| rex_extract | 127 | 121 | 99 |
| rex_extract_all | 153 | 152 | 152 |
| rex_match | 146 | 146 | 124 |
| rex_match_all | 213 | 207 | 183 |

`rex_compile()` takes 0.60–4.56 ms across the eight patterns. Detection on
100 pre-parsed glycans takes 51–105 ms across those patterns.

| Glycans | rex_detect: string input (ms) | Structure + compiled (ms) |
|---:|---:|---:|
| 1 | 3.19 | 0.77 |
| 10 | 17.50 | 7.62 |
| 100 | 91.00 | 75.00 |
| 1,000 | 745.00 | 694.00 |

Precompilation alone saves little in this batch workload. Reusing parsed
structures helps more, especially for small calls. The source shows that
detect/count both collect all matches through `rex_queries()`; extraction
also constructs and canonicalizes subgraphs. These are candidates for
future profiling, not measured attribution of time to individual internals.

## Validation and environment

- Package tests: 259 passed; zero failures, warnings, or skips.
- Each query function's string versus compiled pattern results were checked
  for equality. Counts from string and structure inputs were checked for all
  eight patterns before timing.
- glyrex 0.1.0.9000 at `2a808806c7c88adb9d6eaf7278e5e18ff15ae787`;
  glyrepr 1.1.0.9000; glyparse 0.8.1; igraph 2.3.4; R 4.6.1, macOS arm64.
- `results-2026-10-24/raw.csv`: batch measurements.
- `results-2026-10-24/summary.csv`: median/min/max and per-glycan timing.
- `results-2026-10-24/session.txt`: full session and source revision.

No package functions were modified. Benchmarks are excluded from package builds.
