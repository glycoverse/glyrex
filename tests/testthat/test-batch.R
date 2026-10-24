test_that("duplicate queries preserve recycling, names, and capture schemas", {
  x <- c(
    a = "Gal(b1-4)GlcNAc",
    a = NA_character_,
    b = "Fuc(a1-3)GlcNAc",
    z = "Gal(b1-4)GlcNAc",
    b = "Fuc(a1-3)GlcNAc"
  )
  p <- c(
    "(?<tip>Hex)-HexNAc",
    NA,
    "(?<tip>Fuc)-HexNAc",
    "(?<tip>Hex)-HexNAc",
    "(?<tip>Hex)-HexNAc"
  )
  for (input in list(x, glyrepr::as_glycan_structure(x))) {
    for (f in list(rex_detect, rex_count, rex_extract)) {
      expected <- lapply(seq_along(x), function(i) f(input[i], p[i]))
      if (inherits(expected[[1]], "glyrepr_structure")) {
        expect_equal(
          as.character(f(input, p)),
          stats::setNames(
            vapply(
              expected,
              as.character,
              character(1)
            ),
            names(input)
          )
        )
      } else {
        expect_identical(f(input, p), do.call(c, expected))
      }
    }
    for (f in list(rex_match_all, rex_locate_all, rex_extract_all)) {
      expected <- lapply(seq_along(x), function(i) f(input[i], p[i])[[1L]])
      names(expected) <- names(input)
      if (identical(f, rex_match_all)) {
        expected[[2L]] <- matrix(
          NA_character_,
          1L,
          2L,
          dimnames = list(NULL, c("match", "tip"))
        )
      }
      actual <- f(input, p)
      if (
        identical(f, rex_extract_all) && inherits(input, "glyrepr_structure")
      ) {
        expect_true(all(vapply(
          actual,
          inherits,
          logical(1),
          "glyrepr_structure"
        )))
        actual <- lapply(actual, as.character)
        expected <- lapply(expected, as.character)
      }
      expect_equal(actual, expected)
    }
    expect_identical(rex_which(input, p), which(rex_detect(input, p)))
    expect_equal(rex_subset(input, p), input[c(1L, 3L, 4L)])
  }
  expect_error(
    rex_match(x, c("(?<a>Hex)", "(?<b>Fuc)", NA, "(?<a>Hex)", "(?<b>Fuc)")),
    "identical capture names"
  )
})

test_that("fragment batching preserves substituents, alditols and graph IDs", {
  x <- c(
    "Gal6S(b1-4)[Fuc(a1-3)]GlcNAc-ol",
    "Neu5Ac(a2-3)GalNAc4S6S",
    "Gal(b1-?)GlcNAc",
    "Gal(b1-3/4)GlcNAc"
  )
  g <- glyrepr::as_glycan_structure(x)
  p <- c("Gal6S-(Fuca3-)GlcNAc", "Neu5Ac-GalNAc4S6S", rep("Gal-GlcNAc", 2))
  expect_equal(rex_extract(x, p), as.character(g))
  expect_equal(as.character(rex_extract(g, p)), as.character(g))
  expect_equal(rex_extract(x[1], "Gal6S"), "Gal6S(b1-")
  graph <- as.list(g)[[1]]
  graph <- igraph::permute(graph, rev(seq_len(igraph::vcount(graph))))
  ids <- which(igraph::vertex_attr(graph, "mono") == "Fuc")
  expect_identical(rex_locate(graph, "Fuc")$nodes[[1]], as.integer(ids))
  expect_equal(as.character(rex_extract(graph, "Fuc")), "Fuc(a1-")
})

test_that("query caches do not survive graph or compiled pattern changes", {
  g <- as.list(glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc"))[[1]]
  expect_true(rex_detect(g, "Gal"))
  id <- which(igraph::vertex_attr(g, "mono") == "Gal")
  g <- igraph::set_vertex_attr(g, "mono", index = id, value = "Man")
  expect_false(rex_detect(g, "Gal"))
  expect_true(rex_detect(g, "Man"))
  p <- rex_compile("Gal")
  expect_false(rex_detect(g, p))
  p$tree[[1]]$specs[[1]]$mono <- "Man"
  expect_true(rex_detect(g, p))
  invalid <- igraph::set_vertex_attr(g, "mono", index = id, value = "INVALID")
  expect_error(rex_detect(invalid, "Man"))
})

test_that("fragment serialization preserves ambiguous branch ordering", {
  x <- "Man(b1-?)Gal6S(?1-?)[Man(b1-?)Fuc(?1-?)][Fuc(b1-?)]Man(b1-"
  p <- ".-.-([Hex|Fuc]-)+."
  expected <- c(
    "Man(b1-?)Gal6S(?1-?)[Fuc(?1-?)][Fuc(b1-?)]Man(b1-",
    "Man(b1-?)Fuc(?1-?)[Fuc(b1-?)]Man(b1-"
  )
  expect_identical(rex_extract_all(x, p)[[1]], expected)
  expect_identical(rex_extract(x, p), expected[1])
  expect_identical(unname(rex_match_all(x, p)[[1]][, 1]), expected)
  g <- glyrepr::as_glycan_structure(x)
  expect_identical(unname(rex_match(g, p)[1, 1]), expected[1])
  expect_identical(
    as.character(rex_extract_all(g, p)[[1]]),
    as.character(glyrepr::as_glycan_structure(expected))
  )
})
