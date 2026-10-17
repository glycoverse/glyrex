# Changelog

## glyrex (development version)

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
