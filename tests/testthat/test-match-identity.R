# All permutations also move the root ahead of both competing branches.
identity_permutations <- function(x) {
  if (length(x) == 1L) {
    return(list(x))
  }
  unlist(
    lapply(x, function(i) {
      lapply(identity_permutations(setdiff(x, i)), function(rest) c(i, rest))
    }),
    recursive = FALSE
  )
}

test_that("branch sites survive every node numbering and input branch order", {
  inputs <- c(
    "Gal(b1-2)[GlcNAc(b1-3)][GlcNAc(b1-6)]Man(b1-",
    "GlcNAc(b1-6)[Gal(b1-2)][GlcNAc(b1-3)]Man(b1-"
  )
  pattern <- "Galb2-((?<site>GlcNAc)-)Manb"
  for (x in inputs) {
    expect_equal(rex_count(x, pattern), 2L)
    expect_equal(nrow(rex_match_all(x, pattern)[[1]]), 2L)
    g <- as.list(glyrepr::as_glycan_structure(x))[[1]]
    for (permutation in identity_permutations(seq_len(igraph::vcount(g)))) {
      graph <- igraph::permute(g, permutation)
      ctx <- rex_graph(graph)
      hits <- rex_search(ctx, rex_compile(pattern))
      sites <- which(ctx$mono == "GlcNAc")
      expect_length(hits, 2L)
      expect_setequal(
        vapply(hits, function(hit) hit$captures$site, integer(1)),
        sites
      )
      locations <- rex_locate_all(graph, pattern)[[1]]
      expect_equal(rex_count(graph, pattern), nrow(locations))
      expect_equal(nrow(rex_match_all(graph, pattern)[[1]]), nrow(locations))
      expect_true(rex_detect(graph, pattern))
      expect_equal(rex_locate(graph, pattern), locations[1L, , drop = FALSE])
      for (i in seq_along(hits)) {
        expect_equal(locations$nodes[[i]], sort(hits[[i]]$nodes))
        edges <- igraph::as_edgelist(graph, names = FALSE)
        expect_equal(
          locations$edges[[i]],
          which(
            edges[, 1] %in% hits[[i]]$nodes & edges[, 2] %in% hits[[i]]$nodes
          )
        )
        expect_equal(
          intersect(locations$nodes[[i]], sites),
          hits[[i]]$captures$site
        )
      }
      # The unconstrained path takes the batch fast path in public queries.
      expect_equal(rex_count(graph, "Galb2-(GlcNAc-)Manb"), 2L)
      expect_equal(rex_count(graph, "Galb2-(GlcNAc-)(?<base>Manb)"), 2L)
      expect_equal(rex_count(graph, "GlcNAc-Man"), 2L)
      expect_equal(nrow(rex_locate_all(graph, "GlcNAc-Man")[[1]]), 2L)
      p <- rex_compile("(?<site>GlcNAc)-Man")
      fast <- rex_prepare_pattern(p)
      expect_equal(rex_search(ctx, fast), rex_search(ctx, p))
      expect_equal(
        rex_search_paths(list(ctx, ctx), fast$path),
        rep(list(rex_search(ctx, p)), 2)
      )
      expect_equal(
        rex_search_paths(list(ctx, ctx), fast$path, first = TRUE),
        rep(list(rex_search(ctx, p, first = TRUE)), 2)
      )
    }
  }
})

test_that("identity includes capture assignments, not just the matched nodes", {
  x <- "GlcNAc(b1-3)[GlcNAc(b1-6)]Man(b1-"
  p <- "(?<site>GlcNAc)-(GlcNAc-)Man"
  ctx <- rex_graph(as.list(glyrepr::as_glycan_structure(x))[[1]])
  hits <- rex_search(ctx, rex_compile(p))
  expect_length(hits, 2L)
  expect_equal(sort(hits[[1]]$nodes), sort(hits[[2]]$nodes))
  expect_false(identical(hits[[1]]$captures, hits[[2]]$captures))
  expect_equal(rex_count(x, p), 2L)
  expect_equal(nrow(rex_locate_all(x, p)[[1]]), 2L)
  # Identical serialized fragments/rows still represent two different sites.
  matched <- rex_match_all(x, p)[[1]]
  expect_equal(nrow(matched), 2L)
  expect_equal(matched[1, ], matched[2, ])
  expect_equal(rex_match(x, p), matched[1, , drop = FALSE])
  expect_length(rex_extract_all(x, p)[[1]], 2L)
  expect_equal(rex_count(x, "GlcNAc-(GlcNAc-)Man"), 1L)
  expect_equal(rex_count(x, "[GlcNAc|GlcNAc]-Man"), 2L)
  expect_equal(rex_count("Gal", "(?<site>[Gal|Gal])"), 1L)
})

test_that("subset pruning cannot erase a different capture site", {
  x <- "GlcNAc(b1-4)GlcNAc(b1-4)Man(b1-"
  p <- "(?<site>GlcNAc)-[GlcNAc]?-Man"
  expect_equal(rex_count(x, p), 2L)
  expect_equal(lengths(rex_locate_all(x, p)[[1]]$nodes), c(3L, 2L))
  expect_equal(nrow(rex_match_all(x, p)[[1]]), 2L)
  # With the same capture assignment the nested suffix is still suppressed.
  expect_equal(rex_count(x, "GlcNAc-[GlcNAc]?(?<site>Man)"), 1L)
})

test_that("repeat preference and overlapping matches are retained", {
  x <- "Gal(b1-4)GlcNAc(b1-4)GlcNAc"
  for (suffix in c("*", "*?", "+", "+?")) {
    p <- paste0("Gal-[GlcNAc]", suffix)
    expected <- switch(suffix, "*" = 3L, "*?" = 1L, "+" = 3L, "+?" = 2L)
    expect_equal(rex_count(x, p), 1L)
    expect_equal(lengths(rex_locate_all(x, p)[[1]]$nodes), expected)
    expect_equal(nrow(rex_match_all(x, p)[[1]]), 1L)
  }
  expect_equal(rex_count(x, ".-."), 2L)
  expect_equal(rex_count(x, "[.]{2}"), 2L)
  expect_equal(rex_locate_all(x, ".-."), rex_locate_all(x, "[.]{2}"))
  # Distinct captured lengths are observable identities, even for lazy repeats.
  for (p in c("Gal-(?<site>[GlcNAc]+)", "Gal-(?<site>[GlcNAc]+?)")) {
    expect_equal(rex_count(x, p), 2L)
    expect_equal(lengths(rex_locate_all(x, p)[[1]]$nodes), c(3L, 2L))
  }
})

test_that("identity normalizes node order and capture name order", {
  a <- list(nodes = c(3L, 1L, 2L), captures = list(b = 2L, a = c(3L, 1L)))
  b <- list(nodes = c(1L, 2L, 3L), captures = list(a = c(1L, 3L), b = 2L))
  expect_length(rex_finalize_hits(list(a, b)), 1L)
  b$captures <- list(a = 2L, b = c(1L, 3L))
  expect_length(rex_finalize_hits(list(a, b)), 2L)
  expect_false(identical(
    rex_capture_key(list()),
    rex_capture_key(list(site = integer()))
  ))
})


test_that("batch offsets across decimal boundaries preserve numeric tie order", {
  x <- "Gal(b1-2)[GlcNAc(b1-3)][GlcNAc(b1-6)]Man(b1-"
  g <- as.list(glyrepr::as_glycan_structure(x))[[1]]
  g <- igraph::permute(g, rev(seq_len(igraph::vcount(g))))
  ctx <- rex_graph(g)
  prefix <- rex_graph(as.list(glyrepr::as_glycan_structure(
    "Gal(b1-4)GlcNAc(b1-4)Man"
  ))[[1]])
  contexts <- c(list(prefix), rep(list(ctx), 30L))
  for (pattern in c("GlcNAc-Man", "(?<site>GlcNAc)-Man")) {
    p <- rex_prepare_pattern(rex_compile(pattern))
    for (first in c(FALSE, TRUE)) {
      expect_equal(
        rex_search_paths(contexts, p$path, first = first),
        lapply(contexts, rex_search, compiled = p, first = first)
      )
    }
  }
})
