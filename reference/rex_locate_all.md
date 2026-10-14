# Locate All Structural Matches

Locate the nodes and edges of all matches in each glycan.

## Usage

``` r
rex_locate_all(glycans, pattern)
```

## Arguments

- glycans:

  A character vector of glycan expressions in any format supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html),
  a `glyrepr_structure` vector, or one glyrepr-compatible `igraph`.
  Character inputs may mix formats; malformed structures throw an error.

- pattern:

  A character vector of patterns or one
  [`rex_compile()`](https://glycoverse.github.io/glyrex/reference/rex_compile.md)
  object. Inputs recycle only from length one. Missing glycans or
  patterns propagate.

## Value

A list of data frames with integer list-columns `nodes` and `edges`. No
match produces zero rows; missing input produces one row with
`NA_integer_` in both columns. IDs are one-based and refer to the input
graph, or the canonical graph obtained from character input, as in
[`glyrepr::structure_nodes()`](https://glycoverse.github.io/glyrepr/reference/structure_tables.html)
and
[`glyrepr::structure_edges()`](https://glycoverse.github.io/glyrepr/reference/structure_tables.html).
Only edges within the match are returned.

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

## See also

Other glycan queries:
[`rex_count()`](https://glycoverse.github.io/glyrex/reference/rex_count.md),
[`rex_detect()`](https://glycoverse.github.io/glyrex/reference/rex_detect.md),
[`rex_extract()`](https://glycoverse.github.io/glyrex/reference/rex_extract.md),
[`rex_extract_all()`](https://glycoverse.github.io/glyrex/reference/rex_extract_all.md),
[`rex_locate()`](https://glycoverse.github.io/glyrex/reference/rex_locate.md),
[`rex_match()`](https://glycoverse.github.io/glyrex/reference/rex_match.md),
[`rex_match_all()`](https://glycoverse.github.io/glyrex/reference/rex_match_all.md),
[`rex_subset()`](https://glycoverse.github.io/glyrex/reference/rex_subset.md),
[`rex_which()`](https://glycoverse.github.io/glyrex/reference/rex_which.md)

## Examples

``` r
x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
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
```
