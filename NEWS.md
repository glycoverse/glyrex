# glyrex (development version)

## Breaking changes

* `rex_count()`, `rex_match_all()`, and other queries preserve distinct branch matches and named capture assignments, including captured repeat lengths, regardless of node numbering. Counts can increase and identical extracted strings or location rows can occur; callers needing unique fragments should explicitly deduplicate extracted strings. Greedy/lazy preference and strict-subset filtering now apply only to identical capture assignments, and preference at a starting residue no longer discards incomparable branch choices.

## Minor improvements and fixes

* `rex_match_all()` and `rex_locate_all()` retain distinct valid sites sharing a minimum node ID, with consistent identity and deterministic ordering across general and fixed-path matching.

* `rex_compile()` and query functions are faster for batch queries, including batches of distinct structures (#1).

# glyrex 0.1.0

* Initial release.

* Get Started and glyco-regex syntax vignettes introduce structural queries through biological examples.

* New `rex_compile()` validates and compiles reusable glycan patterns with residue classes, linkage constraints, branches, repetition, anchors, lookaround, and named captures.

* New `rex_count()` and `rex_detect()` count and detect structural matches. Query functions accept glycan strings in any format supported by `glyparse::auto_parse()`, including mixed-format vectors, as well as `glyrepr_structure` vectors and glyrepr-compatible `igraph` objects.

* New `rex_extract()` and `rex_extract_all()` return matching glycan fragments while preserving reducing-end anomeric information.

* New `rex_locate()` and `rex_locate_all()` identify the nodes and edges of structural matches.

* New `rex_match()` and `rex_match_all()` return structural matches and named captures.

* New `rex_subset()` and `rex_which()` select matching glycans or their indices.
