# Subset Glycans by Structural Matches

Keep glycans that match a structural pattern.

## Usage

``` r
rex_subset(glycans, pattern, negate = FALSE)
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

- negate:

  Invert detection when `TRUE`; missing values remain missing.

## Value

A subset preserving the input representation and names. Missing results
are omitted. For a single graph input, returns the graph on a match and
`NULL` otherwise. The pattern must have length one or the length of
`glycans`.

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
[`rex_locate_all()`](https://glycoverse.github.io/glyrex/reference/rex_locate_all.md),
[`rex_match()`](https://glycoverse.github.io/glyrex/reference/rex_match.md),
[`rex_match_all()`](https://glycoverse.github.io/glyrex/reference/rex_match_all.md),
[`rex_which()`](https://glycoverse.github.io/glyrex/reference/rex_which.md)

## Examples

``` r
x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
rex_subset(x, "Gal")
#> [1] "Gal(b1-4)GlcNAc"
rex_subset(x, "Gal", negate = TRUE)
#> [1] "Fuc(a1-3)GlcNAc"
```
