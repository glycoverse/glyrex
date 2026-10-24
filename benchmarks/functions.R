# Run from the package root: Rscript benchmarks/functions.R
# Warm-session elapsed times, including GC; no optional benchmark dependency.
args <- commandArgs(TRUE)
if (length(args) >= 2L) {
  library(glyrex, lib.loc = args[2L])
} else {
  devtools::load_all(quiet = TRUE)
}
set.seed(123)
out <- if (length(args)) {
  args[1L]
} else {
  file.path("benchmarks", "results-current")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fixtures <- dget("tests/testthat/fixtures/glycowork.R")
corpus <- unique(vapply(fixtures, `[[`, character(1), "glycan"))
write.csv(
  data.frame(glycan = corpus),
  file.path(out, "corpus.csv"),
  row.names = FALSE
)
patterns <- c(
  simple = "HexNAc",
  linkage = "Galb4-GlcNAc",
  branch = "Galb4-(Fuca3-)GlcNAc",
  capture = "(?<terminal>Hex)-HexNAc",
  repeat_path = "[HexNAc]{2,}",
  absence = "(!Neu5Ac-)Gal-GlcNAc",
  assertion = "Hex-HexNAc-([Hex|Fuc]-){1,2}HexNAc(?=-HexNAc)",
  no_match = "Kdn"
)
rows <- list()
measure <- function(fun, section, fn, mode, pattern, n) {
  fun()
  loops <- 1L
  repeat {
    t <- system.time(
      for (i in seq_len(loops)) {
        fun()
      }
    )[["elapsed"]]
    if (t >= 0.05 || loops >= 4096L) {
      break
    }
    loops <- loops * 2L
  }
  for (r in seq_len(7L)) {
    gc()
    elapsed <- system.time(
      for (i in seq_len(loops)) {
        fun()
      },
      gcFirst = FALSE
    )[["elapsed"]]
    rows[[length(rows) + 1L]] <<- data.frame(
      section,
      fn,
      mode,
      pattern,
      n,
      replicate = r,
      loops,
      ms = elapsed * 1000 / loops
    )
  }
  cat(section, fn, mode, pattern, n, "done\n")
}
x <- rep(corpus, length.out = 100L)
g <- glyparse::auto_parse(x)
p <- rex_compile(patterns[["capture"]])
fns <- setdiff(getNamespaceExports("glyrex"), "rex_compile")
for (fn in sort(fns)) {
  f <- getExportedValue("glyrex", fn)
  stopifnot(identical(f(x, patterns[["capture"]]), f(x, p)))
  for (mode in c("text", "text_compiled", "structure_compiled")) {
    input <- if (mode == "structure_compiled") g else x
    pat <- if (mode == "text") patterns[["capture"]] else p
    measure(function() f(input, pat), "api", fn, mode, "capture", length(x))
  }
}
for (label in names(patterns)) {
  pat <- patterns[[label]]
  compiled <- rex_compile(pat)
  stopifnot(identical(rex_count(x, pat), rex_count(g, compiled)))
  measure(
    function() rex_compile(pat),
    "compile",
    "rex_compile",
    "text",
    label,
    1L
  )
  measure(
    function() rex_detect(g, compiled),
    "patterns",
    "rex_detect",
    "structure_compiled",
    label,
    length(g)
  )
}
for (n in c(1L, 10L, 100L, 1000L)) {
  x <- rep(corpus, length.out = n)
  g <- glyparse::auto_parse(x)
  for (mode in c("text", "structure_compiled")) {
    input <- if (mode == "text") x else g
    pat <- if (mode == "text") patterns[["capture"]] else p
    measure(
      function() rex_detect(input, pat),
      "scaling",
      "rex_detect",
      mode,
      "capture",
      n
    )
  }
}
# The unique-only control distinguishes engine improvements from deduplication.
x <- corpus
g <- glyparse::auto_parse(x)
for (fn in sort(fns)) {
  f <- getExportedValue("glyrex", fn)
  measure(
    function() f(g, p),
    "unique",
    fn,
    "structure_compiled",
    "capture",
    length(g)
  )
}
raw <- do.call(rbind, rows)
write.csv(raw, file.path(out, "raw.csv"), row.names = FALSE)
keys <- c("section", "fn", "mode", "pattern", "n")
summary <- aggregate(raw["ms"], raw[keys], function(v) {
  c(median = median(v), min = min(v), max = max(v))
})
summary <- cbind(summary[keys], as.data.frame(summary$ms))
summary$ms_per_glycan <- summary$median / summary$n
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
capture.output(
  sessionInfo(),
  Sys.info(),
  if (length(args) >= 3L) args[3L] else system("git rev-parse HEAD", intern = TRUE),
  file = file.path(out, "session.txt")
)
print(summary, row.names = FALSE)
