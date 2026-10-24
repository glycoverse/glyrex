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
        expect_equal(
          rex_search(ctx, fast, first = TRUE),
          utils::head(rex_search(ctx, p), 1L),
          info = pattern
        )
      }
    }
  }
})

test_that("batched paths preserve graph boundaries and original IDs", {
  x <- c(
    "Gal",
    "GlcNAc",
    "Gal(b1-4)GlcNAc",
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc"
  )
  graphs <- as.list(glyrepr::as_glycan_structure(x))
  graphs <- lapply(graphs, function(g) {
    igraph::permute(g, rev(seq_len(igraph::vcount(g))))
  })
  contexts <- lapply(graphs, rex_graph)
  for (pattern in c(
    "Gal-GlcNAc",
    ".",
    ".-.-.",
    "Hex-HexNAc",
    "Kdn",
    "(?<whole>(?<tip>Hex)-(?<base>HexNAc))"
  )) {
    p <- rex_prepare_pattern(rex_compile(pattern))
    expect_identical(
      rex_search_paths(contexts, p$path),
      lapply(contexts, rex_search, compiled = p),
      info = pattern
    )
    expect_identical(
      rex_search_paths(contexts, p$path, first = TRUE),
      lapply(contexts, rex_search, compiled = p, first = TRUE),
      info = pattern
    )
  }
  expect_identical(rex_detect(x[1:2], "Gal-GlcNAc"), c(FALSE, FALSE))
  named <- c(a = "Gal", a = "GlcNAc", b = "Man")
  expect_identical(rex_locate(named, ".")$nodes, I(list(1L, 1L, 1L)))
  expect_identical(
    lapply(rex_locate_all(named, "."), function(x) x$nodes[[1L]]),
    list(a = 1L, a = 1L, b = 1L)
  )
})
