# Compile a Glycan Regular Expression

Compile a reusable structural query. Patterns run from the non-reducing
end towards the reducing end; they are not text regular expressions.

## Usage

``` r
rex_compile(pattern)
```

## Arguments

- pattern:

  One non-missing pattern string, or a compiled pattern.

## Value

A `glyrex_pattern` containing a parsed syntax tree and capture names.

## Details

Residues are separated by `-`. Linkages accept shorthand (`Galb4`) or
IUPAC notation (`Gal(b1-4)`). `.` matches any residue; `Hex`, `HexNAc`,
and other glyrepr classes match their concrete members. `Sia` matches
Neu5Ac, Neu5Gc, and Kdn. `[Gal|Man]` expresses alternatives; `[^Fuc]`
excludes residues. Parentheses describe a sibling branch, as in
`Galb4-(Fuca3)-GlcNAc`. `(!...-)` always asserts branch absence:
`(!Fuc-)` forbids a Fuc branch with any linkage, whereas `(!Fuca3-)`
forbids an alpha1-3 Fuc branch. The assertion checks unconsumed branches
attached to the next residue, consumes no nodes, and can contain a
multi-residue branch pattern. Write the attachment dash inside the
group, immediately before `)`, followed directly by the attachment
residue: `Galb4-(!Fuc-)GlcNAc`. Bare `!residue` and absence groups
without the final dash are rejected.

Bracketed groups accept `?`, `*`, `+`, `{n}`, `{n,m}`, `{n,}`, or
`{,m}`. Append `?` for lazy repetition. Open bounds are limited by the
input graph, not a fixed repeat count. `^` requires a non-reducing
terminus, `$` the reducing end, and `%` an internal residue. Lookahead
`(?=...)` / `(?!...)` tests the reducing-end direction; lookbehind
`(?<=...)` / `(?<!...)` tests the non-reducing direction. Assertions
consume no residues. Empty matches are not returned.

`(?<name>...)` captures a path (an extension to glycowork); `(?:...)`
groups a path without capture. Ordinary parentheses remain branches. The
capture name `match` is reserved for the full-match result column.

An omitted anomer is unconstrained. `GlcNAca` constrains the anomer at
either an internal residue or the reducing end. `GlcNAc(a1-` explicitly
requires the reducing end and its anomeric position, using glyrepr
notation. A closed form, `GlcNAc(a1-)`, is also accepted inside groups.
A specified acceptor position, e.g. `GlcNAca4`, requires an outgoing
glycosidic linkage and cannot match the reducing end. Unknown linkage
fields are wildcards on either side, following glycowork; ambiguous
positions match when their possible values overlap.

## Examples

``` r
p <- rex_compile("Galb4-(Fuca3)-GlcNAc")
rex_detect("Gal(b1-4)[Fuc(a1-3)]GlcNAc", p)
#> [1] TRUE
rex_detect(c("GlcNAc(a1-", "GlcNAc(b1-"), "GlcNAca$")
#> [1]  TRUE FALSE
```
