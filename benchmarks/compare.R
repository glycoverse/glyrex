# Rscript benchmarks/compare.R after running functions.R for both libraries.
before <- read.csv("benchmarks/results-before-optimization/summary.csv")
after <- read.csv("benchmarks/results-after-optimization/summary.csv")
keys <- c("section", "fn", "mode", "pattern", "n")
comparison <- merge(before, after, by = keys, suffixes = c("_before", "_after"))
comparison$speedup <- comparison$median_before / comparison$median_after
write.csv(comparison, "benchmarks/comparison.csv", row.names = FALSE)
table <- function(section, mode = NULL) {
  x <- comparison[comparison$section == section, ]
  if (!is.null(mode)) {
    x <- x[x$mode == mode, ]
  }
  label <- if (section %in% c("compile", "patterns")) "pattern" else "fn"
  if (section == "scaling") {
    x <- x[order(x$n, x$mode), ]
    labels <- paste(x$n, x$mode)
  } else {
    labels <- x[[label]]
  }
  c(
    "| Case | Before (ms) | After (ms) | Speedup |",
    "|---|---:|---:|---:|",
    sprintf(
      "| %s | %.3f | %.3f | %.2fx |",
      labels,
      x$median_before,
      x$median_after,
      x$speedup
    )
  )
}
text <- c(
  "# glyrex performance optimization — 2026-10-24",
  "",
  "Input structure parsing is unchanged: character inputs still call the same `glyparse::auto_parse()` with the same arguments and inputs. No dependency was upgraded. Improvements are in pattern compilation, matching and output construction.",
  "",
  "## Paired measurement",
  "",
  "Old and new packages were installed in separate temporary libraries and measured sequentially in fresh R processes with identical dependencies. Seven timed batches per case; median wall time includes garbage collection during measurement. Package loading, input preparation and precompilation are excluded where indicated. There is no cross-call cache of query results; only immutable residue metadata is cached across calls.",
  "",
  "The shared pattern is `(?<terminal>Hex)-HexNAc`. The 100-input workload repeats 26 fixture glycans (2–15 residues); the unique control uses each of those glycans once. These fixtures are not a population-weighted biological dataset. Structure extraction includes construction of `glyrepr_structure` outputs.",
  "",
  "## 100 structures, precompiled pattern",
  "",
  table("api", "structure_compiled"),
  "",
  "## 100 strings, string pattern (including unchanged parsing)",
  "",
  table("api", "text"),
  "",
  "## 26 distinct structures, precompiled pattern",
  "",
  table("unique"),
  "",
  "## Pattern compilation",
  "",
  table("compile"),
  "",
  "## Detection across pattern types (100 structures)",
  "",
  table("patterns"),
  "",
  "## Detection across batch sizes",
  "",
  table("scaling"),
  "",
  "## Implementation and verification",
  "",
  "- Each unique structure-pattern pair is matched once per call, then results are restored to the original order and names.",
  "- Fixed paths use vectorized matching; general backtracking retains its semantics and reuses atom masks within a query.",
  "- Raw igraph inputs retain full validation. Existing structure vectors reuse their validated graphs; floating metadata remains rejected.",
  "- Fragments are serialized in batches; ambiguous branch-order ties use the original graph extraction path to preserve exact strings.",
  "- 790 baseline-versus-optimized result objects agree exactly after normalizing structure objects to their typed, named IUPAC values. They cover all 10 query functions, 39 patterns, 26 fixture glycans, 30 seeded random glycans (2–8 residues), missing values, and permuted raw graphs.",
  "- Package tests: 485 passed. R CMD check: 0 errors, 0 warnings, 0 notes (`--no-manual`, `_R_CHECK_SYSTEM_CLOCK_=FALSE` to disable unavailable remote clock verification).",
  "- Fixed-path tests include 192 assertions comparing the fast path with the backtracking engine. A discovered ambiguous-branch serialization difference has a dedicated regression test.",
  "- All 39 compiled syntax trees are identical to the baseline after the final token extraction change.",
  "",
  "## Reproduce",
  "",
  "Install each source revision into a separate R library, then run:",
  "",
  "```sh",
  "Rscript benchmarks/functions.R benchmarks/results-before-optimization /path/to/before-library BEFORE_REVISION",
  "Rscript benchmarks/functions.R benchmarks/results-after-optimization /path/to/after-library AFTER_REVISION",
  "Rscript benchmarks/compare.R",
  "Rscript benchmarks/parity.R /path/to/before-library /tmp/before.rds",
  "Rscript benchmarks/parity.R /path/to/after-library /tmp/after.rds",
  "Rscript -e 'stopifnot(identical(readRDS(\"/tmp/before.rds\"), readRDS(\"/tmp/after.rds\")))'",
  "```",
  "",
  "Raw batch timings, corpus copies, dependency versions and source revisions are retained in `results-before-optimization/` and `results-after-optimization/`. `comparison.csv` contains all 64 paired cases. The original first-run baseline remains in `results-2026-10-24/`."
)
writeLines(text, "benchmarks/optimization.md")
print(
  comparison[
    comparison$section == "api",
    c("fn", "mode", "median_before", "median_after", "speedup")
  ],
  row.names = FALSE
)
