# glyrex 0.0.0.9000

* `(!...-)` always asserts branch absence: `(!Fuc-)` forbids any Fuc branch at the attachment residue, while `(!Fuca3-)` restricts the forbidden branch to an alpha1-3 linkage.

* `rex_compile()` compiles glycowork-style glycan patterns with branches, repetition, location anchors, lookaround, and named captures, including glyrepr reducing-end configuration constraints.
* `rex_detect()`, `rex_count()`, `rex_extract()`, `rex_extract_all()`, `rex_match()`, `rex_match_all()`, `rex_locate()`, `rex_locate_all()`, `rex_subset()`, and `rex_which()` provide stringr-style structural queries on IUPAC strings and glyrepr structures.
* Queries reject glycans with floating parts or floating substituents.
