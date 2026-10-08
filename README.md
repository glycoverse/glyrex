
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

`glyrex` is being developed as a unified structural query language for
glycans within the [glycoverse](https://github.com/glycoverse)
ecosystem. It combines a [stringr](https://stringr.tidyverse.org/)-style
API, a glycan-aware regular expression language, and a graph-native
matching engine specialized for rooted glycan trees.

Where stringr queries strings with regular expressions, glyrex aims to
query glycan structures with glyco-regular expressions. Patterns
describe residues, glycosidic linkages, paths, and branches rather than
positions in an IUPAC string.

**Development status:** The package is currently a scaffold. The API,
pattern syntax, and examples below describe the intended design and are
not implemented or runnable yet.

## Installation

Currently, glyrex is available as a development version from
[GitHub](https://github.com/glycoverse/glyrex):

``` r
# install.packages("pak")
pak::pkg_install("glycoverse/glyrex")
```

CRAN, r-universe, and GitHub release installation instructions will be
added when those distributions become available.

**Note:** General installation tips for the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse) are also useful
here: [Installation of
glycoverse](https://github.com/glycoverse/glycoverse#installation).

## Documentation

Function reference documentation and a getting-started guide will be
added as the public API develops. For now, follow development and report
issues on [GitHub](https://github.com/glycoverse/glyrex).

## Role in `glycoverse`

[glyrepr](https://github.com/glycoverse/glyrepr) provides glycan
representations; `glyrex` will provide the language and operations for
searching those structures. Character inputs will be converted through
glyrepr, while existing glyrepr structure objects will be matched
directly.

The planned query objects will complement the motif analysis tools in
[glymotif](https://github.com/glycoverse/glymotif) and provide a
reusable structural query interface for other glycoverse packages.

## Planned API

The API follows familiar stringr operations, adapted to glycan
structures:

| Function                              | Intended operation                                             |
|:--------------------------------------|:---------------------------------------------------------------|
| `rex_detect()`                        | Test whether each glycan contains a pattern.                   |
| `rex_count()`                         | Count structural matches in each glycan.                       |
| `rex_extract()` / `rex_extract_all()` | Extract the first or all matching subglycans.                  |
| `rex_match()` / `rex_match_all()`     | Return matches and named captures.                             |
| `rex_locate()` / `rex_locate_all()`   | Locate matching nodes and edges.                               |
| `rex_subset()`                        | Keep glycans that match a pattern.                             |
| `rex_which()`                         | Return the indices of matching glycans.                        |
| `rex_compile()`                       | Compile a pattern once for reuse as a `glyrex_pattern` object. |
| `rex_validate()` / `rex_explain()`    | Validate patterns and explain their parsed structure.          |

Extraction is intended to preserve the input representation: character
input returns character output, and glyrepr input returns glyrepr
structures. Locations will identify graph nodes and edges rather than
character offsets. Pattern objects (`glyrex_pattern`) and match objects
(`glyrex_match`) will be separate, so queries can be reused and results
can retain structural context.

## Planned pattern language

Patterns will support residue names and shared residue classes, linkage
constraints, wildcards, sequences, branches, alternatives, negation,
quantifiers, anchors, and lookaround. Named captures are planned
separately from ordinary parentheses, which express branches.

For example, `Neu5Aca3-Galb4-GlcNAc` is intended to describe a
structural path with specified residues and linkages. Generic residue
classes such as `Hex` and `HexNAc` will use shared glycoverse residue
definitions.

Patterns will be parsed into an explicit abstract syntax tree and
compiled for matching against rooted glycan trees. Open-ended
quantifiers (`*`, `+`, and `{n,}`) will be bounded by the input
structure, with no fixed repeat limit.

The initial scope targets rooted, connected glycan trees. Nested branch
patterns, repeat units, ambiguous topology, composition-only inputs, and
graph rewriting are deferred. Capture support is planned beyond the
initial core search operations.

## Example

The following illustrates the planned workflow. These functions are not
yet available, and this example is not executed when building the
README.

``` r
library(glyrex)

glycans <- c(
  "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
  "Gal(b1-4)GlcNAc"
)

# Compile once and reuse across operations
pattern <- rex_compile("Neu5Aca3-Galb4-GlcNAc")

rex_detect(glycans, pattern)
rex_subset(glycans, pattern)
rex_count(glycans, pattern)
rex_extract_all(glycans, pattern)

# Locate matches by glycan node and edge IDs
rex_locate_all(glycans, pattern)
```

A future named-capture query could use
`(?<antenna>Neu5Aca3-Galb4)-GlcNAc` with `rex_match_all()` to extract
both the full match and the captured antenna.
