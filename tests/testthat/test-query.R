test_that("vectorized APIs preserve missing, duplicate, and empty inputs", {
  x <- c(a = "Gal(b1-4)GlcNAc", b = "Fuc(a1-3)GlcNAc", a = NA_character_)
  expect_equal(rex_detect(x, "Gal"), c(a = TRUE, b = FALSE, a = NA))
  expect_equal(rex_count(x, "Gal"), c(a = 1L, b = 0L, a = NA_integer_))
  expect_equal(rex_subset(x, "Gal"), x[1])
  expect_equal(rex_subset(x, "Gal", negate = TRUE), x[2])
  expect_equal(unname(rex_which(x, "Gal")), 1L)
  expect_equal(
    rex_extract_all(x, "Gal"),
    list(a = "Gal(b1-", b = character(), a = NA_character_)
  )
  expect_equal(rex_detect(character(), "Gal"), logical())
  expect_equal(rex_count("Gal", character()), integer())
  expect_equal(rex_extract(character(), "Gal"), character())
  expect_equal(nrow(rex_locate(character(), "Gal")), 0L)
  expect_equal(rex_detect("Gal", c("Gal", "Man", NA)), c(TRUE, FALSE, NA))
  expect_equal(
    rex_detect(x, c("Gal", "Fuc", "Man")),
    c(a = TRUE, b = TRUE, a = NA)
  )
  expect_equal(rex_locate_all(x, "Gal")[[2]]$nodes, I(list()))
  expect_equal(rex_locate_all(x, "Gal")[[3]]$nodes, I(list(NA_integer_)))
  expect_equal(rex_detect(x, rex_compile("Gal")), rex_detect(x, "Gal"))
})

test_that("captures return matrices and preserve absent group columns", {
  x <- c("Gal(b1-4)GlcNAc(a1-", "Fuc(a1-3)GlcNAc", NA)
  p <- "(?<tip>Gal)-(?<root>GlcNAc)"
  m <- rex_match(x, p)
  expect_equal(dim(m), c(3L, 3L))
  expect_equal(colnames(m), c("match", "tip", "root"))
  expect_equal(unname(m[1, ]), c(x[1], "Gal(b1-", "GlcNAc(a1-"))
  expect_equal(unname(m[2:3, ]), matrix(NA_character_, 2, 3))
  all <- rex_match_all(x, p)
  expect_equal(vapply(all, nrow, integer(1)), c(1L, 0L, 1L))
  expect_equal(
    unname(rex_match("GlcNAc", "(?<tip>Gal)?-GlcNAc")[1, "tip"]),
    NA_character_
  )
  expect_equal(
    unname(rex_match("Gal(b1-4)GlcNAc", "(?<both>(?<tip>Gal)-GlcNAc)")[
      1,
      "tip"
    ]),
    "Gal(b1-"
  )
  s <- glyrepr::as_glycan_structure(x)
  expect_s3_class(rex_extract(s, "Gal"), "glyrepr_structure")
  expect_s3_class(rex_extract_all(s, "Gal")[[2]], "glyrepr_structure")
  expect_identical(rex_subset(s, "Gal"), s[1])
})

test_that("unsupported inputs and incompatible recycling fail clearly", {
  expect_snapshot(
    error = TRUE,
    rex_detect(c("Gal", "Man"), c("Gal", "Man", "Fuc"))
  )
  expect_snapshot(error = TRUE, rex_subset("Gal", c("Gal", "Man")))
  expect_snapshot(error = TRUE, rex_detect(1, "Gal"))
  expect_snapshot(error = TRUE, rex_detect("Gal", "Gal", negate = NA))
  expect_snapshot(error = TRUE, rex_detect("{Fuc(a1-3)}Gal(b1-4)GlcNAc", "Gal"))
})


test_that("empty matching retains the compiled capture schema", {
  expect_equal(dim(rex_match(character(), "(?<tip>Gal)")), c(0L, 2L))
  expect_equal(
    colnames(rex_match(character(), "(?<tip>Gal)")),
    c("match", "tip")
  )
})

test_that("floating substituents are rejected for every representation", {
  x <- glyrepr::as_glycan_structure("{6S|1,2}Gal(b1-4)GlcNAc")
  expect_snapshot(error = TRUE, rex_detect(x, "Gal"))
  expect_snapshot(error = TRUE, rex_count(as.list(x)[[1]], "Gal"))
  expect_snapshot(error = TRUE, rex_extract(as.character(x), "Gal"))
})
