# Console-only comparison; neither outputs nor datasets are saved.
source("benchmarks/helpers.R")
before <- load_glyrex_version("2a808806c7c88adb9d6eaf7278e5e18ff15ae787")
after <- load_glyrex_version()
set.seed(123)
fixtures <- dget("tests/testthat/fixtures/glycowork.R")
generated <- lapply(seq_len(30), function(i) {
  n <- sample(2:8, 1)
  parent <- vapply(2:n, function(j) sample.int(j - 1L, 1), integer(1))
  # Unknown acceptor positions allow arbitrary branches without site conflicts.
  mono <- sample(c("Gal", "Man", "Fuc", "GlcNAc", "GalNAc", "Neu5Ac"), n, TRUE)
  list(
    mono = mono,
    sub = sample(c("", "", "", "6S"), n, TRUE),
    edges = as.integer(rbind(parent, 2:n)),
    linkage = paste0(
      sample(c("a", "b", "?"), n - 1L, TRUE),
      ifelse(mono[-1L] == "Neu5Ac", "2", "1"),
      "-?"
    ),
    anomer = paste0(
      sample(c("a", "b", "?"), 1),
      ifelse(mono[1] == "Neu5Ac", "2", "1")
    ),
    alditol = i %% 5L == 0L
  )
})
x <- c(
  unique(vapply(fixtures, `[[`, character(1), "glycan")),
  as.character(glyrepr::structure_from_arrays(generated)),
  NA_character_
)
names(x) <- rep(c("a", "a", "b", ""), length.out = length(x))
patterns <- unique(c(
  vapply(fixtures, `[[`, character(1), "pattern"),
  "(?<tip>Hex)-HexNAc",
  "(?<arm>(?<tip>Gal)-GlcNAc)",
  "[.]+",
  ".",
  "(?<tip>Gal)?-GlcNAc",
  "Gal6S",
  "[HexNAc]{2,}",
  "(!Neu5Ac-)Gal-GlcNAc"
))
compared <- 0L
g <- glyrepr::as_glycan_structure(x)
selected <- commandArgs(TRUE)
if (length(selected)) query_functions <- intersect(query_functions, selected)
for (fn in query_functions) {
  for (j in seq_along(patterns)) {
    for (representation in c("text", "structure")) {
      input <- if (representation == "text") x else g
      key <- paste(fn, j, representation, sep = "/")
      p <- before$rex_compile(patterns[j])
      stopifnot(identical(p, after$rex_compile(patterns[j])))
      old <- normalize_result(before[[fn]](input, p))
      new <- normalize_result(after[[fn]](input, p))
      if (!identical(old, new)) {
        cat("Mismatch:", key, "pattern:", patterns[j], "\n")
        print(all.equal(old, new))
        stop("Parity check failed.")
      }
      compared <- compared + 1L
    }
  }
  # Raw graphs retain original node and edge IDs after a permutation.
  graph <- igraph::permute(
    as.list(g)[[1]],
    rev(seq_len(igraph::vcount(as.list(g)[[1]])))
  )
  stopifnot(identical(
    normalize_result(before[[fn]](graph, "(?<tip>Hex)-HexNAc")),
    normalize_result(after[[fn]](graph, "(?<tip>Hex)-HexNAc"))
  ))
  compared <- compared + 1L
  cat(fn, "complete\n")
}
cat(compared, "result objects agree\n")
