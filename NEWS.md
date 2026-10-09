# glyrex 0.0.0.9000

* Query functions now use `glycans` instead of `string` and automatically parse all formats supported by `glyparse::auto_parse()`, including mixed-format vectors.

* Branch absence requires `(!branch-)`; bare `!linked-residue` and groups missing the final attachment dash are rejected. `[^A]` always consumes a nonmatching residue, including when A has a linkage constraint.

* Residue negation uses `[^A]`; the `[!A]` spelling is rejected with migration guidance.

* `(!...-)` always asserts branch absence: `(!Fuc-)` forbids any Fuc branch at the attachment residue, while `(!Fuca3-)` restricts the forbidden branch to an alpha1-3 linkage.

* `rex_compile()` compiles glycowork-style glycan patterns with branches, repetition, location anchors, lookaround, and named captures, including glyrepr reducing-end configuration constraints.
* `rex_detect()`, `rex_count()`, `rex_extract()`, `rex_extract_all()`, `rex_match()`, `rex_match_all()`, `rex_locate()`, `rex_locate_all()`, `rex_subset()`, and `rex_which()` provide stringr-style structural queries on IUPAC strings and glyrepr structures.
* Queries reject glycans with floating parts or floating substituents.
