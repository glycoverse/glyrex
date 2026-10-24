# Changelog

## glyrex (development version)

### Breaking changes

- [`rex_count()`](https://glycoverse.github.io/glyrex/reference/rex_count.md),
  [`rex_match_all()`](https://glycoverse.github.io/glyrex/reference/rex_match_all.md),
  and other queries preserve distinct branch matches and named capture
  assignments, including captured repeat lengths, regardless of node
  numbering. Counts can increase and identical extracted strings or
  location rows can occur; callers needing unique fragments should
  explicitly deduplicate extracted strings. Greedy/lazy preference and
  strict-subset filtering now apply only to identical capture
  assignments, and preference at a starting residue no longer discards
  incomparable branch choices.

### Minor improvements and fixes

- [`rex_match_all()`](https://glycoverse.github.io/glyrex/reference/rex_match_all.md)
  and
  [`rex_locate_all()`](https://glycoverse.github.io/glyrex/reference/rex_locate_all.md)
  retain distinct valid sites sharing a minimum node ID, with consistent
  identity and deterministic ordering across general and fixed-path
  matching.

- [`rex_compile()`](https://glycoverse.github.io/glyrex/reference/rex_compile.md)
  and query functions are faster for batch queries, including batches of
  distinct structures (#1).

## glyrex 0.1.0

- Initial release.

- Get Started and glyco-regex syntax vignettes introduce structural
  queries through biological examples.

- New
  [`rex_compile()`](https://glycoverse.github.io/glyrex/reference/rex_compile.md)
  validates and compiles reusable glycan patterns with residue classes,
  linkage constraints, branches, repetition, anchors, lookaround, and
  named captures.

- New
  [`rex_count()`](https://glycoverse.github.io/glyrex/reference/rex_count.md)
  and
  [`rex_detect()`](https://glycoverse.github.io/glyrex/reference/rex_detect.md)
  count and detect structural matches. Query functions accept glycan
  strings in any format supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html),
  including mixed-format vectors, as well as `glyrepr_structure` vectors
  and glyrepr-compatible `igraph` objects.

- New
  [`rex_extract()`](https://glycoverse.github.io/glyrex/reference/rex_extract.md)
  and
  [`rex_extract_all()`](https://glycoverse.github.io/glyrex/reference/rex_extract_all.md)
  return matching glycan fragments while preserving reducing-end
  anomeric information.

- New
  [`rex_locate()`](https://glycoverse.github.io/glyrex/reference/rex_locate.md)
  and
  [`rex_locate_all()`](https://glycoverse.github.io/glyrex/reference/rex_locate_all.md)
  identify the nodes and edges of structural matches.

- New
  [`rex_match()`](https://glycoverse.github.io/glyrex/reference/rex_match.md)
  and
  [`rex_match_all()`](https://glycoverse.github.io/glyrex/reference/rex_match_all.md)
  return structural matches and named captures.

- New
  [`rex_subset()`](https://glycoverse.github.io/glyrex/reference/rex_subset.md)
  and
  [`rex_which()`](https://glycoverse.github.io/glyrex/reference/rex_which.md)
  select matching glycans or their indices.
