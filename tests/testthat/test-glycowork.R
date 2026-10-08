# Cases adapted from BojarLab/glycowork tests/test_core_functions.py at
# 8e178129e99477dfdc38a055a34886bbbe0e02ea (MIT; see fixtures/LICENSE-glycowork).
# Excludes compact-IUPAC input conversion and explicitly unsupported floating
# inputs. Compare canonical structures after neutralizing only root anomers:
# glycowork omits them, while glyrex deliberately preserves them.
test_that("glycowork reference cases retain structural semantics", {
  cases <- dget(test_path("fixtures", "glycowork.R"))
  normalize <- function(x) {
    if (!length(x)) {
      return(character())
    }
    graphs <- as.list(glyrepr::as_glycan_structure(x))
    sort(vapply(
      graphs,
      function(g) {
        g <- igraph::set_graph_attr(g, "anomer", "?1")
        glyrepr::graph_to_iupac(glyrepr::canonicalize_glycan_graph(g))
      },
      character(1)
    ))
  }
  for (case in cases) {
    actual <- rex_extract_all(case$glycan, case$pattern)[[1]]
    info <- paste("Upstream line", case$line, "pattern", case$pattern)
    if (is.logical(case$expected)) {
      expect_equal(length(actual) > 0L, case$expected, info = info)
    } else {
      expect_equal(
        normalize(actual),
        normalize(unlist(case$expected)),
        info = info
      )
    }
    expect_equal(
      rex_detect(case$glycan, case$pattern),
      length(actual) > 0L,
      info = info
    )
  }
})

test_that("lookbehinds constrain individual branches", {
  p <- "(?<!HexNAc-)Mana3-(?<!HexNAc-)([Mana6]){1}-Manb4-GlcNAcb4-GlcNAc"
  x <- c(
    "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    "Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
  )
  expect_equal(rex_detect(x, p), c(TRUE, FALSE, FALSE))
})

test_that("negated linked groups assert absence without consuming residues", {
  p <- "Galb3-(!GlcNAcb6-)GalNAc"
  expect_equal(
    rex_detect(c("Gal(b1-3)GalNAc", "Gal(b1-3)[GlcNAc(b1-6)]GalNAc"), p),
    c(TRUE, FALSE)
  )
})

test_that("optional sialic acid and sibling branches compose with lookahead", {
  p <- "r[Sia]{,1}-Monosaccharide-([dHex]){,1}-Monosaccharide(?=-Mana6-Monosaccharide)"
  x <- "GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Gc(a2-6)GalNAc(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
  expected <- glyrepr::as_glycan_structure(
    "Neu5Gc(a2-6)GalNAc(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(rex_extract(x, p), as.character(expected))
})
