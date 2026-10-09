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

test_that("queries accept named glycans in every glyparse format", {
  glycans <- c(
    condensed = "Gal(b1-4)GlcNAc(b1-",
    extended = "beta-D-Galp-(1->3)-alpha-D-GalpNAc-(1->",
    glycoct = paste(
      "RES",
      "1b:a-dgal-HEX-1:5",
      "2s:n-acetyl",
      "3b:b-dgal-HEX-1:5",
      "LIN",
      "1:1d(2+1)2n",
      "2:1o(3+1)3d",
      sep = "\n"
    ),
    short = "Neu5Aca3Gala3(Fuca6)GlcNAcb-",
    glycam = "DManpa1-3[DManpa1-6]DManpb1-4DGlcpNAcb1-4DGlcpNAcb1-OH",
    compact = "Mana1-3(Mana1-6)Manb1-4GlcNAcb",
    wurcs = paste0(
      "WURCS=2.0/2,3,2/",
      "[a2112h-1b_1-5][a2112h-1a_1-4]/1-2-2/a4-b1_b2-c1"
    ),
    linear_code = "Ma3(Ma6)Mb4GNb4GNb",
    pglyco = "(N(F)(N(H(H(N))(H(N(H))))))",
    strucgp = "A2B2C1D1E2F1fedD1E2edcbB5ba",
    kcf = paste(
      "ENTRY       G00001                      Glycan",
      "NODE        2",
      "            1   Glc        0     0",
      "            2   Gal        6     0",
      "EDGE        1",
      "            1     2:b1    1:4",
      "///",
      sep = "\n"
    ),
    linucs = "[][b-D-Glcp]{[(4+1)][b-D-Galp]{}}",
    gwb = paste0(
      "freeEnd--1b1D-GlcNAc,p(--6a1L-Fuc,p)",
      "--4b1D-Gal,p--3a2D-NeuAc,p$MONO,Und,0,0,freeEnd"
    )
  )
  generic <- paste0(
    "Hex(??-?)HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]",
    "Hex(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-"
  )
  canonical <- c(
    "Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-3)GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(a1-3)[Fuc(a1-6)]GlcNAc(b1-",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-",
    "Galf(a1-2)Galf(a1-4)Gal(b1-",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    generic,
    generic,
    "Gal(b1-4)Glc(?1-",
    "Gal(b1-4)Glc(b1-",
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
  )
  names(canonical) <- names(glycans)
  glycans <- c(glycans, missing = NA_character_, glycans[1])
  canonical <- c(canonical, missing = NA_character_, canonical[1])
  queries <- list(
    rex_detect,
    rex_count,
    rex_extract,
    rex_extract_all,
    rex_match,
    rex_match_all,
    rex_locate,
    rex_locate_all
  )
  for (query in queries) {
    expect_identical(
      query(glycans = glycans, pattern = "."),
      query(glycans = canonical, pattern = ".")
    )
  }
  expect_identical(
    rex_count(glycans = glycans, pattern = "."),
    setNames(
      c(2L, 2L, 2L, 4L, 5L, 4L, 3L, 5L, 9L, 9L, 2L, 2L, 4L, NA_integer_, 2L),
      names(glycans)
    )
  )
  expect_identical(
    rex_subset(glycans = glycans, pattern = "."),
    glycans[!is.na(glycans)]
  )
  expect_identical(
    rex_which(glycans = glycans, pattern = "."),
    which(!is.na(glycans))
  )
  expect_identical(
    rex_detect(glycans = c(NA_character_, NA_character_), "."),
    c(NA, NA)
  )
})

test_that("malformed expressions fail through the automatic parser", {
  expect_snapshot(
    error = TRUE,
    rex_detect(glycans = c("Galb1-4GlcNAc", "WURCS-invalid"), pattern = "Gal")
  )
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
