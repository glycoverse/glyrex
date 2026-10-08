# Query Glycan Structures with Regular Expressions

A stringr-like interface to graph-native glycan pattern matching.

## Usage

``` r
rex_detect(string, pattern, negate = FALSE)

rex_count(string, pattern)

rex_extract(string, pattern)

rex_extract_all(string, pattern)

rex_match(string, pattern)

rex_match_all(string, pattern)

rex_locate(string, pattern)

rex_locate_all(string, pattern)

rex_which(string, pattern, negate = FALSE)

rex_subset(string, pattern, negate = FALSE)
```

## Arguments

- string:

  An IUPAC-condensed character vector, a `glyrepr_structure` vector, or
  one glyrepr-compatible `igraph`. Character inputs are parsed by
  glyrepr; malformed structures throw an error.

- pattern:

  A character vector of patterns or one
  [`rex_compile()`](https://glycoverse.github.io/glyrex/reference/rex_compile.md)
  object. Inputs recycle only from length one. Missing strings or
  patterns propagate.

- negate:

  Invert detection when `TRUE`; missing values remain missing.

## Value

`rex_detect()` returns logicals; `rex_count()` returns integers.
`rex_extract()` returns the first match, or `NA`; `rex_extract_all()`
returns a list of vectors (empty for no match, length-one `NA` for
missing). Extraction returns canonical glyrepr IUPAC strings for
character input and glyrepr structure vectors for structure or graph
input. A cut fragment retains the anomer and donor position of its
former attachment linkage.

`rex_match()` returns a character matrix with `match` followed by named
capture columns; unmatched captures are `NA`. `rex_match_all()` returns
a list of such matrices, with zero rows for no match and one missing row
for missing input. Vectorized patterns must have identical capture
names.

`rex_locate()` returns a data frame with integer list-columns `nodes`
and `edges`; `rex_locate_all()` returns a list of these data frames. IDs
are one-based and refer to the input graph, or the canonical graph
obtained from character input, as in
[`glyrepr::structure_nodes()`](https://glycoverse.github.io/glyrepr/reference/structure_tables.html)
and
[`glyrepr::structure_edges()`](https://glycoverse.github.io/glyrepr/reference/structure_tables.html).
Only edges within the match are returned. First locations use
`NA_integer_` for no match or missing input; all locations use zero rows
for no match and one missing row for missing input.

`rex_subset()` preserves the input representation and names;
`rex_which()` returns matching indices. Both omit missing results and
require `pattern` to have length one or the length of `string`.

## Details

The matching unit is a connected structural trace. As in glycowork, the
first successful trace per starting node is retained, then strict
subsets of other traces are removed. Distinct branches with equal
sequences remain distinct matches. Matches are ordered by decreasing
residue count and then starting node ID. Repeats are greedy unless
marked lazy; successful traces may overlap. Counts therefore need not
equal text-regex counts. Only connected rooted trees are supported;
floating topology or floating substituents are rejected explicitly. See
[`rex_compile()`](https://glycoverse.github.io/glyrex/reference/rex_compile.md)
for syntax.

## Examples

``` r
x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
rex_detect(x, "Gal-HexNAc")
#> [1]  TRUE FALSE    NA
rex_count(x, "HexNAc")
#> [1]  1  1 NA
rex_extract_all(x, "HexNAc")
#> [[1]]
#> [1] "GlcNAc(?1-"
#> 
#> [[2]]
#> [1] "GlcNAc(?1-"
#> 
#> [[3]]
#> [1] NA
#> 
rex_match(x, "(?<terminal>Gal)-GlcNAc")
#>      match                 terminal 
#> [1,] "Gal(b1-4)GlcNAc(?1-" "Gal(b1-"
#> [2,] NA                    NA       
#> [3,] NA                    NA       
rex_locate_all(x, "Gal-GlcNAc")
#> [[1]]
#>   nodes edges
#> 1  1, 2     1
#> 
#> [[2]]
#> [1] nodes edges
#> <0 rows> (or 0-length row.names)
#> 
#> [[3]]
#>   nodes edges
#> 1    NA    NA
#> 
rex_subset(x, "Gal")
#> [1] "Gal(b1-4)GlcNAc"
rex_which(x, "Fuc")
#> [1] 2
```
