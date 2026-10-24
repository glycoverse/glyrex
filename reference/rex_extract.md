# Extract the First Structural Match

Extract the first matching fragment from each glycan.

## Usage

``` r
rex_extract(glycans, pattern)
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

A character vector of canonical glyrepr IUPAC strings for character
input, or a glyrepr structure vector for structure or graph input. No
match or missing input returns a missing value. A cut fragment retains
the anomer and donor position of its former attachment linkage.

## Details

A match is identified by its complete set of graph nodes and the node
set assigned to each named capture. Repeated search routes to the same
identity are collapsed; different capture assignments remain distinct,
even when their extracted strings or location rows are identical.

At the same actual starting residue (the first consumed node), the first
successful trace wins over later traces with identical capture
assignments whose node sets contain it or are contained in it.
Incomparable branch choices are retained. Search order follows pattern
alternatives and greedy (longer first) or lazy (shorter first)
repetition. Among these selected traces, strict subsets are removed only
when capture assignments agree. Thus different captured sites or
captured repeat lengths are preserved, including for lazy repeats. An
unmatched capture differs from an empty one.

Results are ordered by decreasing residue count, then minimum node ID,
then lexicographically by sorted node IDs and capture assignments to
break ties. This order is deterministic for a given input graph;
renumbering can change the order, but IDs always refer to that input
graph. Matches may overlap. Counts therefore need not equal text-regex
counts. These rules replace the earlier one-match-per-minimum-node and
capture-blind subset filtering. Only connected rooted trees are
supported; floating topology or floating substituents are rejected
explicitly. See
[`rex_compile()`](https://glycoverse.github.io/glyrex/reference/rex_compile.md)
for syntax.

## See also

Other glycan queries:
[`rex_count()`](https://glycoverse.github.io/glyrex/reference/rex_count.md),
[`rex_detect()`](https://glycoverse.github.io/glyrex/reference/rex_detect.md),
[`rex_extract_all()`](https://glycoverse.github.io/glyrex/reference/rex_extract_all.md),
[`rex_locate()`](https://glycoverse.github.io/glyrex/reference/rex_locate.md),
[`rex_locate_all()`](https://glycoverse.github.io/glyrex/reference/rex_locate_all.md),
[`rex_match()`](https://glycoverse.github.io/glyrex/reference/rex_match.md),
[`rex_match_all()`](https://glycoverse.github.io/glyrex/reference/rex_match_all.md),
[`rex_subset()`](https://glycoverse.github.io/glyrex/reference/rex_subset.md),
[`rex_which()`](https://glycoverse.github.io/glyrex/reference/rex_which.md)

## Examples

``` r
x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
rex_extract(x, "HexNAc")
#> [1] "GlcNAc(?1-" "GlcNAc(?1-" NA          
```
