# glyrex

`glyrex` searches glycan structures using a glycowork-style regular
expression language and a stringr-style R API. It provides a more
flexible and expressive way for glycan substructure analysis than
[glymotif](https://github.com/glycoverse/glymotif).

## Installation

The glyrex package is not a core glycoverse package. You need to install
it individually even if you have installed the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse).

You can install the latest release of glyrex from
[r-universe](https://glycoverse.r-universe.dev/glyrex)
(**recommended**):

``` r

pak::repo_add(glycoverse = "https://glycoverse.r-universe.dev")
pak::pkg_install("glyrex")
```

Or from [GitHub](https://github.com/glycoverse/glyrex):

Install the latest releases of this package and its required Glycoverse
dependencies together from GitHub:

``` r

pak::pkg_install(c(
  "glycoverse/glyrepr@*release",
  "glycoverse/glyparse@*release",
  "glycoverse/glyrex@*release"
))
```

Or install the development version (NOT recommended):

``` r

pak::pkg_install("glycoverse/glyrex")
```

**Note:** Tips and troubleshooting for the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse) are also
applicable here: [Installation of
glycoverse](https://github.com/glycoverse/glycoverse#installation).

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
