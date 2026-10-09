
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glyrex <a href="https://glycoverse.github.io/glyrex/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-universe
version](https://glycoverse.r-universe.dev/glyrex/badges/version)](https://glycoverse.r-universe.dev/glyrex)
[![R-CMD-check](https://github.com/glycoverse/glyrex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/glycoverse/glyrex/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/glycoverse/glyrex/graph/badge.svg)](https://app.codecov.io/gh/glycoverse/glyrex)
<!-- badges: end -->

`glyrex` searches glycan structures using a glycowork-style regular
expression language and a stringr-style R API. Matching follows glycan
graph topology, including branches, rather than character positions in
an IUPAC string.

## Installation

``` r
# install.packages("pak")
pak::pkg_install("glycoverse/glyrex")
```

This is an experimental development package. It depends on
[glyrepr](https://glycoverse.github.io/glyrepr/) for glycan
representations.

## Role in glycoverse

`glymotif` provides a concise syntax for matching glycan motifs.
However, exact motif matching may not be sufficient to express more
complex structural constraints. For example, a sialyl-LacNAc motif
without fucosylation cannot be specified using a bare IUPAC-condensed
sequence. To address this limitation, glyrex introduces a dedicated
glyco-regex syntax that enables more flexible and expressive pattern
matching.

## Example

``` r
library(glyrex)
x <- c(
  "Gal(b1-4)[Fuc(a1-3)]GlcNAc(a1-",
  "Gal(b1-4)GlcNAc(b1-",
  NA_character_
)
p <- rex_compile("Galb4-(Fuca3-)GlcNAc")
rex_detect(x, p)
#> [1]  TRUE FALSE    NA
rex_extract(x, p)
#> [1] "Fuc(a1-3)[Gal(b1-4)]GlcNAc(a1-" NA
#> [3] NA
rex_match(x, "(?<antenna>Galb4)-GlcNAc")
#>      match                 antenna
#> [1,] "Gal(b1-4)GlcNAc(a1-" "Gal(b1-"
#> [2,] "Gal(b1-4)GlcNAc(b1-" "Gal(b1-"
#> [3,] NA                    NA
rex_locate_all(x, p)
#> [[1]]
#>     nodes edges
#> 1 1, 2, 3  1, 2
#>
#> [[2]]
#> [1] nodes edges
#> <0 rows> (or 0-length row.names)
#>
#> [[3]]
#>   nodes edges
#> 1    NA    NA
```

## Acknowledgement

The glyco-regex syntax is adapted from
[glycowork](https://github.com/BojarLab/glycowork/), while the API
design is inspired by [stringr](https://github.com/tidyverse/stringr).
