# glyrex (development version)

* `rex_compile()` and query functions are faster for batch queries, including batches of distinct structures.

# glyrex 0.1.0

* Initial release.

* Get Started and glyco-regex syntax vignettes introduce structural queries through biological examples.

* New `rex_compile()` validates and compiles reusable glycan patterns with residue classes, linkage constraints, branches, repetition, anchors, lookaround, and named captures.

* New `rex_count()` and `rex_detect()` count and detect structural matches. Query functions accept glycan strings in any format supported by `glyparse::auto_parse()`, including mixed-format vectors, as well as `glyrepr_structure` vectors and glyrepr-compatible `igraph` objects.

* New `rex_extract()` and `rex_extract_all()` return matching glycan fragments while preserving reducing-end anomeric information.

* New `rex_locate()` and `rex_locate_all()` identify the nodes and edges of structural matches.

* New `rex_match()` and `rex_match_all()` return structural matches and named captures.

* New `rex_subset()` and `rex_which()` select matching glycans or their indices.
