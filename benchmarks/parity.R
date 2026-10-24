# Run in separate processes with old/new library paths and compare the RDS files.
# Rscript benchmarks/parity.R /path/to/library /path/to/output.rds
args <- commandArgs(TRUE)
library(glyrex, lib.loc = args[1L])
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
normalize <- function(value) {
  if (inherits(value, "glyrepr_structure")) {
    return(list(
      type = "glyrepr_structure",
      values = as.character(value)
    ))
  }
  if (inherits(value, "igraph")) {
    return(glyrepr::structure_to_iupac(value))
  }
  if (is.list(value) && !is.data.frame(value)) {
    return(lapply(value, normalize))
  }
  value
}
out <- list()
g <- glyrepr::as_glycan_structure(x)
for (fn in setdiff(getNamespaceExports("glyrex"), "rex_compile")) {
  f <- getExportedValue("glyrex", fn)
  for (j in seq_along(patterns)) {
    for (representation in c("text", "structure")) {
      input <- if (representation == "text") x else g
      key <- paste(fn, j, representation, sep = "/")
      out[[key]] <- normalize(f(input, rex_compile(patterns[j])))
    }
  }
  # Raw graphs retain original node and edge IDs after a permutation.
  graph <- igraph::permute(
    as.list(g)[[1]],
    rev(seq_len(igraph::vcount(as.list(g)[[1]])))
  )
  out[[paste0(fn, "/graph")]] <- normalize(f(graph, "(?<tip>Hex)-HexNAc"))
  cat(fn, "complete\n")
}
saveRDS(out, args[2L])
cat(length(out), "result objects saved\n")
