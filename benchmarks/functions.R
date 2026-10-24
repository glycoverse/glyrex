# Run from the package root: Rscript benchmarks/functions.R [baseline-ref]
# Timings and comparisons live only in memory and are printed to the console.
source("benchmarks/helpers.R")
args <- commandArgs(TRUE)
baseline_ref <- if (length(args)) {
  args[1L]
} else {
  "2a808806c7c88adb9d6eaf7278e5e18ff15ae787"
}
before <- load_glyrex_version(baseline_ref)
after <- load_glyrex_version()
fixtures <- dget("tests/testthat/fixtures/glycowork.R")
x <- unique(vapply(fixtures, `[[`, character(1), "glycan"))
g <- glyparse::auto_parse(x)
stopifnot(!anyDuplicated(as.character(g)))
p <- before$rex_compile("(?<terminal>Hex)-HexNAc")
measure <- function(f) {
  for (i in 1:5) {
    f(g, p)
  }
  loops <- 1L
  repeat {
    elapsed <- system.time(
      for (i in seq_len(loops)) {
        f(g, p)
      }
    )[["elapsed"]]
    if (elapsed >= 0.1) {
      break
    }
    loops <- loops * 2L
  }
  vapply(
    seq_len(9L),
    function(i) {
      gc()
      system.time(
        for (j in seq_len(loops)) {
          f(g, p)
        },
        gcFirst = FALSE
      )[["elapsed"]] *
        1000 /
        loops
    },
    numeric(1)
  )
}
rows <- lapply(seq_along(query_functions), function(i) {
  fn <- query_functions[i]
  stopifnot(identical(
    normalize_result(before[[fn]](g, p)),
    normalize_result(after[[fn]](g, p))
  ))
  # Alternate which version is timed first across functions.
  if (i %% 2L) {
    old <- measure(before[[fn]])
    new <- measure(after[[fn]])
  } else {
    new <- measure(after[[fn]])
    old <- measure(before[[fn]])
  }
  data.frame(
    fn,
    before_ms = median(old),
    after_ms = median(new),
    speedup = median(old) / median(new)
  )
})
cat("Unique structures:", length(g), "\nBaseline:", baseline_ref, "\n")
print(do.call(rbind, rows), row.names = FALSE)
