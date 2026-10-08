test_that("paths and branches match graph topology", {
  x <- "Gal(b1-4)[Fuc(a1-3)]GlcNAc(a1-"
  expect_equal(
    rex_detect(
      x,
      c(
        "Galb4-(Fuca3)-GlcNAc",
        "Fuca3-(Galb4)-GlcNAc",
        "Gal-Fuc",
        "Galb3-GlcNAc"
      )
    ),
    c(TRUE, TRUE, FALSE, FALSE)
  )
  expect_equal(rex_count(x, "."), 3L)
  expect_equal(rex_count(x, "[Gal|Fuc]"), 2L)
  expect_equal(rex_count(x, "Hex-HexNAc"), 1L)
  expect_equal(rex_count(x, "[^Fuc]"), 2L)
  expect_equal(rex_locate(x, "Galb4-(Fuca3)-GlcNAc")$nodes[[1]], 1:3)
  expect_equal(length(rex_locate(x, "Galb4-(Fuca3)-GlcNAc")$edges[[1]]), 2L)
})

test_that("reducing-end anomers and extracted fragment anomers are preserved", {
  x <- c("Gal(b1-4)GlcNAc(a1-", "Gal(b1-4)GlcNAc(b1-", "Gal(b1-4)GlcNAc(?1-")
  expect_equal(rex_detect(x, "GlcNAc$"), rep(TRUE, 3))
  expect_equal(rex_detect(x, "GlcNAca$"), c(TRUE, FALSE, TRUE))
  expect_equal(rex_detect(x, "GlcNAc(a1-"), c(TRUE, FALSE, TRUE))
  expect_equal(rex_detect(x, "GlcNAc(b1-)"), c(FALSE, TRUE, TRUE))
  expect_equal(rex_detect(x, "GlcNAca4$"), rep(FALSE, 3))
  expect_equal(rex_extract(x, "Gal"), rep("Gal(b1-", 3))
  expect_equal(
    rex_extract(x, "GlcNAc"),
    c("GlcNAc(a1-", "GlcNAc(b1-", "GlcNAc(?1-")
  )
  expect_equal(
    rex_detect(
      "Neu5Ac(a2-3)Gal",
      c("Neu5Aca3-Gal", "Neu5Ac(a1-3)-Gal", "Sia-Gal")
    ),
    c(TRUE, FALSE, TRUE)
  )
  expect_equal(rex_extract("Gal(b1-4)GlcNAc-ol", "GlcNAc"), "GlcNAc-ol(?1-")
})

test_that("repetition is bounded by the graph and respects greediness", {
  x <- "Gal(b1-4)GlcNAc(b1-4)GlcNAc"
  expect_equal(rex_extract(x, "Hex-[HexNAc]*?"), "Gal(b1-")
  expect_equal(rex_extract(x, "Hex-[HexNAc]+?"), "Gal(b1-4)GlcNAc(b1-")
  expect_equal(rex_detect(x, "[HexNAc]{3}"), FALSE)
  expect_equal(rex_detect(x, "[HexNAc]{2,}"), TRUE)
  expect_equal(rex_detect(x, "[Hex]{,2}"), TRUE)
  expect_equal(rex_detect("GlcNAc", "[Gal]?"), FALSE)
  long <- paste0(paste(rep("Gal(b1-4)", 12), collapse = ""), "GlcNAc")
  expect_equal(length(rex_locate(long, "[Gal]+-GlcNAc")$nodes[[1]]), 13L)
  expect_equal(rex_detect("Gal", "[Gal]{999999999}"), FALSE)
})

test_that("assertions consume no nodes and bind to the right direction", {
  x <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(a1-"
  expect_equal(
    rex_detect(x, c("^Gal", "^Neu5Ac", "Gal%", "GlcNAc%", "Gal$")),
    c(FALSE, TRUE, TRUE, FALSE, FALSE)
  )
  expect_equal(rex_extract(x, "(?<=Neu5Ac-)Gal(?=-GlcNAc)"), "Gal(b1-")
  expect_equal(rex_detect(x, "(?<!Neu5Ac-)Gal"), FALSE)
  expect_equal(rex_detect(x, "Gal(?!-Man)"), TRUE)
  expect_equal(rex_detect(x, "Gal(?!-GlcNAc)"), FALSE)
  expect_equal(rex_count(x, "(?=Gal)"), 0L)
  expect_equal(
    rex_detect(
      c("Gal(b1-4)GlcNAc", "Gal(b1-4)[Fuc(a1-3)]GlcNAc"),
      "Gal-!Fuca3-GlcNAc"
    ),
    c(TRUE, FALSE)
  )
})

test_that("linkage ambiguity and residue modifications follow structural matching", {
  expect_equal(rex_detect("Gal(b1-?)GlcNAc", "Galb3/4-GlcNAc"), TRUE)
  expect_equal(rex_detect("Gal(b1-3/4)GlcNAc", "Galb4-GlcNAc"), TRUE)
  expect_equal(rex_detect("Gal(b1-3)GlcNAc", "Galb4-GlcNAc"), FALSE)
  expect_equal(
    rex_detect("Gal6S(b1-4)GlcNAc", c("Gal6S", "Gal", ".")),
    c(TRUE, FALSE, TRUE)
  )
  expect_equal(rex_detect("D-IdoA(b1-3)GalNAc4S", "D-IdoA-GalNAc4S"), TRUE)
})

test_that("locations reference the supplied graph without renumbering", {
  g <- as.list(glyrepr::as_glycan_structure("Gal(b1-4)[Fuc(a1-3)]GlcNAc"))[[1]]
  g <- igraph::permute(g, c(3L, 1L, 2L))
  id <- which(igraph::vertex_attr(g, "mono") == "Gal")
  expect_equal(rex_locate(g, "Gal")$nodes[[1]], id)
  expect_equal(rex_locate(g, "Gal")$edges[[1]], integer())
  expect_s3_class(rex_extract(g, "Gal"), "glyrepr_structure")
})


test_that("branch-only tails cannot yield disconnected full matches", {
  x <- "Gal(b1-4)[Fuc(a1-3)]GlcNAc"
  expect_equal(rex_detect(x, "Gal-(Fuc)"), FALSE)
  expect_equal(rex_extract_all(x, "Gal-(Fuc)")[[1]], character())
  expect_equal(rex_detect(x, "Gal-(Fuc)-GlcNAc"), TRUE)
})

test_that("parenthesized negation always asserts absence at the attachment", {
  x <- c(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-6)]GlcNAc",
    "Neu5Ac(a2-3)Gal(b1-4)[Man(a1-6)]GlcNAc",
    "Neu5Ac(a2-3)[Fuc(a1-2)]Gal(b1-4)GlcNAc",
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(b1-3)]GlcNAc"
  )
  any_fuc <- "Neu5Aca3-Galb4-(!Fuc)-GlcNAc"
  a3_fuc <- "Neu5Aca3-Galb4-(!Fuca3)-GlcNAc"
  expect_equal(rex_detect(x, any_fuc), c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(rex_detect(x, a3_fuc), c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(rex_detect(x, rex_compile(any_fuc)), rex_detect(x, any_fuc))
  expect_equal(rex_count(x, any_fuc), c(1L, 0L, 0L, 1L, 1L, 0L))
  expect_equal(
    rex_extract(x[4], any_fuc),
    rex_extract(x[4], "Neu5Aca3-Galb4-GlcNAc")
  )
  expect_equal(
    rex_locate(x[4], any_fuc),
    rex_locate(x[4], "Neu5Aca3-Galb4-GlcNAc")
  )
  expect_equal(rex_detect("Fuc(a1-3)GlcNAc", "Fuc-(!Fuc)-GlcNAc"), TRUE)
})

test_that("absence groups also accept full linkages and branch paths", {
  x <- c("Gal(b1-4)GlcNAc", "Gal(b1-4)[Fuc(a1-3)]GlcNAc")
  expect_equal(rex_detect(x, "Gal-(!Fuc(a1-3))-GlcNAc"), c(TRUE, FALSE))
  expect_equal(rex_detect(x, "(!Fuc)-GlcNAc"), c(TRUE, FALSE))
  expect_equal(
    rex_detect("Gal(b1-4)[Fuc(a1-2)Man(a1-6)]GlcNAc", "Gal-(!Fuc-Man)-GlcNAc"),
    FALSE
  )
  expect_equal(
    rex_detect("Gal(b1-4)[Man(a1-6)]GlcNAc", "Gal-(!Fuc-Man)-GlcNAc"),
    TRUE
  )
  expect_equal(rex_detect("GlcNAc", "(!Fuc)"), FALSE)
})
