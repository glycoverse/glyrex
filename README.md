
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

## Example

``` r
library(glyrex)
x <- c(
  "Gal(b1-4)[Fuc(a1-3)]GlcNAc(a1-",
  "Gal(b1-4)GlcNAc(b1-",
  NA_character_
)
p <- rex_compile("Galb4-(Fuca3)-GlcNAc")
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

## API

| stringr             | glyrex              | Result                                           |
|:--------------------|:--------------------|:-------------------------------------------------|
| `str_detect()`      | `rex_detect()`      | Logical vector.                                  |
| `str_count()`       | `rex_count()`       | Number of structural matches.                    |
| `str_extract()`     | `rex_extract()`     | First matching structure.                        |
| `str_extract_all()` | `rex_extract_all()` | List of all matching structures.                 |
| `str_match()`       | `rex_match()`       | Character matrix: full match and named captures. |
| `str_match_all()`   | `rex_match_all()`   | List of match matrices.                          |
| `str_locate()`      | `rex_locate()`      | First matching node and edge IDs.                |
| `str_locate_all()`  | `rex_locate_all()`  | List of all matching locations.                  |
| `str_subset()`      | `rex_subset()`      | Matching input glycans.                          |
| `str_which()`       | `rex_which()`       | Matching input indices.                          |
| —                   | `rex_compile()`     | Reusable `glyrex_pattern` object.                |

Character inputs are IUPAC-condensed strings. `glyrepr_structure`
vectors and single glyrepr-compatible `igraph` objects are also
accepted. Extraction preserves the representation: character input
returns canonical IUPAC strings; structure or graph input returns
glyrepr structure vectors. Match matrices always contain canonical IUPAC
strings. Subsetting preserves the original input.

Patterns can be vectors, with length-one recycling. Incompatible lengths
are errors. `rex_subset()` and `rex_which()` require one pattern or one
per input glycan. Missing input or patterns propagate as `NA`; filtering
omits them. All-match functions distinguish no match (empty result) from
missing input (a length-one missing result). Invalid structures are
errors.

Locations contain `nodes` and `edges` integer list-columns. IDs match
`glyrepr::structure_nodes()` and `glyrepr::structure_edges()` on the
input structure, or on the canonical glyrepr conversion of character
input. Direct graph inputs retain their existing vertex and edge IDs.

## Pattern language

| Syntax                 | Meaning                                                 | Example                   |
|:-----------------------|:--------------------------------------------------------|:--------------------------|
| Residue / class        | A residue or shared glyrepr residue class               | `Gal`, `HexNAc`, `Sia`    |
| `.`                    | Any residue                                             | `.-GlcNAc`                |
| `-`                    | Connected path towards the reducing end                 | `Gal-GlcNAc`              |
| Linkage                | Shorthand or full linkage constraint                    | `Galb4`, `Gal(b1-4)`      |
| `[...]`                | Alternatives separated by a vertical bar                | See below.                |
| `[^A]`, `[!A]`         | Negated residue                                         | `[^Fuc]-GlcNAc`           |
| `(path)`               | Sibling branch at the next attachment residue           | `Galb4-(Fuca3)-GlcNAc`    |
| `!linked-residue`      | Assert absence of an unconsumed branch                  | `Galb3-!GlcNAcb6-GalNAc`  |
| `?`, `*`, `+`, `{n,m}` | Repeat bracketed units or branches                      | `[HexNAc]{2,}`            |
| Lazy suffix `?`        | Prefer fewer repetitions                                | `Hex-[HexNAc]+?`          |
| `^`, `$`, `%`          | Non-reducing terminus, reducing end, internal position  | `^Gal`, `GlcNAc$`, `Hex%` |
| `(?=...)`, `(?!...)`   | Positive / negative lookahead towards reducing end      | `Gal(?=-GlcNAc)`          |
| `(?<=...)`, `(?<!...)` | Positive / negative lookbehind towards non-reducing end | `(?<=Sia-)Gal`            |
| `(?<name>...)`         | Named capture of a connected path                       | `(?<arm>Gal-GlcNAc)-Man`  |
| `(?:...)`              | Non-capturing path group                                | `(?:Gal-GlcNAc)+`         |

For example, either galactose or mannose linked to GlcNAc:

``` r
rex_compile("[Gal|Man]-GlcNAc")
```

Ordinary parentheses retain the glycowork branch meaning; they do not
create captures. Capture names must be unique (`match` is reserved), and
vectorized match patterns must have the same capture names. An omitted
or optional capture returns `NA`. The legacy leading `r` pattern marker
is accepted.

Unknown linkage fields are compatible with specified fields on either
side, following glycowork. Slash-separated positions match when their
possible values overlap. Substituted residues are matched explicitly
(`Gal6S`); `Gal` does not match `Gal6S`. Classes follow glyrepr, with
the glycowork alias `Sia` for `Neu5Ac`, `Neu5Gc`, and `Kdn`.

Repetition is bounded by the size of the input graph. The engine keeps
the first successful trace per starting node, removes traces strictly
contained in larger ones, and orders results by decreasing residue count
and then starting node ID. Distinct branches can yield identical
extracted sequences. Matches can overlap; counts are structural, not
text-regex occurrence counts. Zero-length matches are omitted.

## Reducing-end configuration

Unlike glycowork, glyrepr records the reducing-end anomer and donor
position. An unqualified residue leaves these unconstrained. A suffix
such as `GlcNAca` constrains the anomer of either an internal residue or
the reducing end. `GlcNAc(a1-` explicitly requires an alpha reducing end
with donor position 1; inside groups, use the closed spelling
`GlcNAc(a1-)`. An acceptor constraint such as `GlcNAca4` requires a real
linkage and cannot match the reducing end.

``` r
ends <- c("GlcNAc(a1-", "GlcNAc(b1-")
rex_detect(ends, "GlcNAc$")
#> [1] TRUE TRUE
rex_detect(ends, "GlcNAca$")
#> [1]  TRUE FALSE
rex_detect(ends, "GlcNAc(a1-")
#> [1]  TRUE FALSE
rex_extract("Gal(b1-4)GlcNAc(a1-", "Gal")
#> [1] "Gal(b1-"
```

An extracted internal fragment retains its former attachment’s anomer
and donor position at its new reducing end. A fragment containing the
original reducing end retains that end’s configuration and alditol
status.

## Scope and compatibility

Only connected rooted glycan trees are accepted. **Glycans containing
floating parts or floating substituents are rejected**, for character,
structure, and graph inputs. Composition-only inputs and repeat-unit
topology are unsupported.

The semantic reference is the [glycowork regex
module](https://github.com/BojarLab/glycowork/blob/8e178129e99477dfdc38a055a34886bbbe0e02ea/glycowork/motif/regex.py).
The test suite includes 49 applicable cases adapted from its upstream
tests, plus branch assertions, named captures, reducing-end
configuration, vector behavior, and floating-input rejection. Upstream
output comparisons normalize only reducing-end configuration and
canonical branch order. Compact-IUPAC input conversion and floating
structures are outside the supported input contract. Unknown residue
names are compilation errors. Named captures and non-capturing path
groups extend the upstream language; open repetition bounds use the
graph size instead of a fixed limit.
