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
