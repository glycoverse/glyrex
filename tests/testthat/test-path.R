test_that("fixed paths agree with backtracking on canonical and permuted graphs", {
  x <- c(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc",
    "Gal6S(b1-3/4)[Fuc(a1-3)]GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-?)GlcNAc(a1-",
    "Hex(??-?)HexNAc"
  )
  patterns <- c(
    ".",
    ".-.",
    ".-.-.",
    "Hex-HexNAc",
    "Galb4-GlcNAc",
    "Galb3/4-GlcNAc",
    "Siaa3-Gal",
    "Gal6S",
    "Monosaccharide",
    "(?<whole>(?<tip>Hex)-(?<base>HexNAc))",
    "[Gal]-GlcNAc",
    "Kdn"
  )
  for (graph in as.list(glyrepr::as_glycan_structure(x))) {
    for (g in list(
      graph,
      igraph::permute(graph, rev(seq_len(igraph::vcount(graph))))
    )) {
      ctx <- rex_graph(g)
      for (pattern in patterns) {
        p <- rex_compile(pattern)
        fast <- rex_prepare_pattern(p)
        expect_false(is.null(fast$path))
        expect_equal(rex_search(ctx, fast), rex_search(ctx, p), info = pattern)
      }
    }
  }
})
