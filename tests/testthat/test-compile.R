test_that("compilation is reusable and captures have explicit syntax", {
  p <- rex_compile("(?<arm>Galb4)-([Fuca3]){1}-GlcNAc")
  expect_s3_class(p, "glyrex_pattern")
  expect_identical(rex_compile(p), p)
  expect_equal(p$captures, "arm")
  expect_output(print(p), "<glyrex_pattern>")
  expect_equal(rex_compile("(Fuc)-GlcNAc")$captures, character())
})

test_that("malformed patterns fail at compilation", {
  expect_snapshot(error = TRUE, rex_compile("[Gal"))
  expect_snapshot(error = TRUE, rex_compile("[Gal]{3,2}"))
  expect_snapshot(error = TRUE, rex_compile("[Gal|]"))
  expect_snapshot(error = TRUE, rex_compile("Gal*"))
  expect_snapshot(error = TRUE, rex_compile("(?<x>Gal)-(?<x>Man)"))
  expect_snapshot(error = TRUE, rex_compile("[GlcNAc]{,}"))
  expect_snapshot(error = TRUE, rex_compile("(?=)"))
  expect_snapshot(error = TRUE, rex_compile(NA_character_))
})


test_that("the full-match column name cannot be shadowed by a capture", {
  expect_snapshot(error = TRUE, rex_compile("(?<match>Gal)"))
})

test_that("bracketed exclamation negation is rejected", {
  expect_snapshot(error = TRUE, rex_compile("[!Fuc]"))
  expect_snapshot(error = TRUE, rex_compile("Gal-([!Fuca3])-GlcNAc"))
  expect_snapshot(error = TRUE, rex_compile("[Gal|!Man]"))
  expect_equal(rex_detect(c("Fuc", "Gal"), "[^Fuc]"), c(FALSE, TRUE))
  expect_equal(rex_detect("Gal(b1-4)GlcNAc", "Gal-(!Fuc-)GlcNAc"), TRUE)
})

test_that("only explicit branch absence groups are accepted", {
  expect_snapshot(error = TRUE, rex_compile("Gal-!Fuca3-GlcNAc"))
  expect_snapshot(error = TRUE, rex_compile("!Fuc"))
  expect_snapshot(error = TRUE, rex_compile("Gal-(!Fuc)-GlcNAc"))
  expect_snapshot(error = TRUE, rex_compile("Gal-(!Fuca3)-GlcNAc"))
  expect_equal(rex_detect("Gal(b1-4)GlcNAc", "Gal-(!Fuca3-)GlcNAc"), TRUE)
  expect_equal(rex_detect("Gal(b1-4)GlcNAc", "[^Fuca3]-GlcNAc"), TRUE)
  expect_equal(rex_detect("GlcNAc", "[^Fuca3]-GlcNAc"), FALSE)
  expect_equal(rex_detect("Fuc(a1-3)GlcNAc", "[^Fuca3]-GlcNAc"), FALSE)
})

test_that("patterns reject IUPAC-style linkages", {
  expect_snapshot(error = TRUE, rex_compile("Gal(b1-4)GlcNAc"))
  expect_snapshot(error = TRUE, rex_compile("Gal(b1-4)-GlcNAc"))
  expect_snapshot(error = TRUE, rex_compile("GlcNAc(a1-"))
  expect_snapshot(error = TRUE, rex_compile("[GlcNAc(b1-)]"))
  expect_snapshot(error = TRUE, rex_compile("Gal-(!Fuc(a1-3)-)GlcNAc"))
  expect_snapshot(error = TRUE, rex_compile("[Gal(?1-3/4)]-GlcNAc"))
  expect_snapshot(
    error = TRUE,
    rex_detect("Gal(b1-4)GlcNAc", "Gal(b1-4)GlcNAc")
  )
  expect_equal(rex_detect("Gal(b1-4)GlcNAc", "Galb4-GlcNAc"), TRUE)
})
