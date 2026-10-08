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
