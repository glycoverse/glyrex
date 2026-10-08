# Changelog

## glyrex 0.0.0.9000

- Branch absence requires `(!branch-)`; bare `!linked-residue` and
  groups missing the final attachment dash are rejected. `[^A]` always
  consumes a nonmatching residue, including when A has a linkage
  constraint.

- Residue negation uses `[^A]`; the `[!A]` spelling is rejected with
  migration guidance.

- `(!...-)` always asserts branch absence: `(!Fuc-)` forbids any Fuc
  branch at the attachment residue, while `(!Fuca3-)` restricts the
  forbidden branch to an alpha1-3 linkage.

- [`rex_compile()`](https://glycoverse.github.io/glyrex/reference/rex_compile.md)
  compiles glycowork-style glycan patterns with branches, repetition,
  location anchors, lookaround, and named captures, including glyrepr
  reducing-end configuration constraints.

- [`rex_detect()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  [`rex_count()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  [`rex_extract()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  [`rex_extract_all()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  [`rex_match()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  [`rex_match_all()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  [`rex_locate()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  [`rex_locate_all()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  [`rex_subset()`](https://glycoverse.github.io/glyrex/reference/rex_query.md),
  and
  [`rex_which()`](https://glycoverse.github.io/glyrex/reference/rex_query.md)
  provide stringr-style structural queries on IUPAC strings and glyrepr
  structures.

- Queries reject glycans with floating parts or floating substituents.
