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

### Reading a structural pattern

A pattern describes connected residues in a glycan graph. Read a path
from left to right, towards the reducing end. Separate adjacent residues
with `-`: for example, `Gal-GlcNAc` asks for a Gal attached directly to
a GlcNAc. It does not skip intervening residues. Without end anchors,
this path can occur within a larger glycan and does not have to describe
the whole structure.

Start with residue names, then add constraints as needed:

- `Gal` matches a galactose residue.

- `.` matches any one residue, not any number of residues.

- `Hex` and `HexNAc` match members of the corresponding glyrepr classes;
  for example, `Hex` matches Gal, Glc, and Man.

- `Sia` matches Neu5Ac, Neu5Gc, and Kdn.

- `[Gal|Man]` matches either Gal or Man. Alternatives can also be paths,
  as in `[Gal-GlcNAc|Gal-GalNAc]`.

- `[^Fuc]` matches one residue that is not Fuc. Negated alternatives
  must each describe a single residue; `[^Fuc|Xyl]` excludes both Fuc
  and Xyl.

### Specifying linkages and the reducing end

A linkage belongs to the residue on its non-reducing side. For example,
`Galb4-GlcNAc` requires Gal to be attached to GlcNAc by a beta1-4
linkage. Use glyco-regex linkage suffixes such as `Galb4`; IUPAC-style
linkages such as `Gal(b1-4)` are not supported in patterns. The suffix
specifies the anomer and acceptor position. The donor position is 1, or
2 for Neu5Ac, Neu5Gc, Kdn, and `Sia`.

Omitted linkage fields impose no constraint. Thus `Gal-GlcNAc` does not
specify a linkage, and `GlcNAca` specifies only the alpha anomer, either
at an internal residue or at the reducing end. A specified acceptor
position, as in `GlcNAca4`, requires an outgoing glycosidic linkage
towards the reducing end and cannot match the reducing-end residue
itself.

To require reducing-end GlcNAc with an alpha anomer, use `GlcNAca$`.

Unknown linkage fields (`?`) act as wildcards in either the pattern or
the input glycan, following glycowork. A specified linkage can therefore
match an input whose corresponding field is unknown; a match establishes
compatibility, not certainty about that linkage. Ambiguous positions
match when their possible values overlap, so acceptor `3/4` is
compatible with `4` but not `6`.

### Requiring or excluding branches

Ordinary parentheses describe a sibling branch attached to the next
residue in the main path. The attachment dash must be inside the
parentheses, with no dash after the branch group. In
`Galb4-(Fuca3-)GlcNAc`, both Gal and Fuc attach to the same GlcNAc: Gal
by beta1-4 and Fuc by alpha1-3. The Fuc is not inserted between Gal and
GlcNAc. Required branch residues are part of the match. Additional
branches are allowed unless explicitly excluded. Pattern parentheses
serve a different purpose from the square brackets used for branches in
an input IUPAC string such as `Gal(b1-4)[Fuc(a1-3)]GlcNAc`.

Use `(!...-)` to require the absence of a sibling branch. For example,
`Galb4-(!Fuc-)GlcNAc` matches the Gal-GlcNAc path only if the GlcNAc has
no remaining Fuc branch with any linkage. `Galb4-(!Fuca3-)GlcNAc`
excludes only a branch compatible with alpha1-3 Fuc. The check considers
branches not already used by the preceding pattern and can describe
several residues, as in `Gal-(!Gal-GlcNAc-)Man`.

Keep the attachment dash inside the absence group, immediately before
`)`, and write the attachment residue directly after it. Bare `!Fuc`,
`[!Fuc]`, and absence groups without the final dash are rejected.

Unlike `[^Fuc]`, which matches and includes one non-Fuc residue,
`(!Fuc-)` only checks a condition and adds no residues to the match.
Consequently, `[^Fuc]-GlcNAc` requires a non-Fuc residue before GlcNAc;
`(!Fuc-)GlcNAc` can match GlcNAc alone if it has no Fuc branch.

### Grouping and repetition

Use square brackets or `(?:...)` to group a path. Ordinary parentheses
mean a branch, not a path group. Attach a repetition operator to the
group:

- `[Gal]?`: zero or one Gal.

- `[Gal]*`: zero or more consecutive Gal residues.

- `[Gal]+`: one or more consecutive Gal residues.

- `[Gal]{2}`: exactly two consecutive Gal residues.

- `[Gal]{2,4}`: between two and four consecutive Gal residues.

- `[Gal]{2,}`: at least two consecutive Gal residues.

- `[Gal]{,4}`: zero to four consecutive Gal residues.

A group can contain more than one residue: `[Galb4-GlcNAcb3]{2}` repeats
the two-residue path twice. Use grouped syntax for repetition rather
than forms such as `Gal*`, which are rejected.

Repetition is greedy by default: longer repetitions are tried first.
Append `?` to try shorter repetitions first, as in `[Gal]+?`. In either
case, the rest of the pattern must also match. Open upper bounds are
limited by the input graph, not a fixed repeat count. A pattern that
permits zero repetitions never produces an empty match.

### Anchors and lookaround assertions

Assertions restrict where a match is allowed without adding residues:

- `^Gal` requires Gal at a non-reducing terminus (a residue with no
  further residues attached on its non-reducing side).

- `GlcNAc$` requires GlcNAc at the reducing end.

- `Gal%` requires Gal to be internal: neither a non-reducing terminus
  nor the reducing-end residue.

Lookaround checks a neighbouring path without including it in the
result. Lookahead examines the reducing-end direction: `Gal-(?=GlcNAc)`
matches Gal only when the next residue is GlcNAc, whereas
`Gal-(?!GlcNAc)` requires that the next residue does not match GlcNAc
(or is absent). Lookbehind examines the non-reducing direction:
`(?<=Gal-)GlcNAc` matches GlcNAc when a Gal is attached to it, and
`(?<!Gal-)GlcNAc` requires that no such Gal path exists. Paths inside
lookbehind are still written towards the reducing end, just like the
rest of the pattern.

Negative lookbehind and branch absence are equivalent at the start of a
pattern: `(?<!Gal-)GlcNAc` and `(!Gal-)GlcNAc` both match only GlcNAc
with no directly attached Gal on its non-reducing side, regardless of
linkage. Neither assertion adds residues to the match.

After a preceding path, the two assertions differ. Negative lookbehind
checks all attached paths, including residues already matched, whereas
branch absence checks only paths not already used by the preceding
pattern. For example, on `Gal(b1-4)GlcNAc`, `Gal-(!Gal-)GlcNAc` matches
because there is no additional Gal branch. `Gal-(?<!Gal-)GlcNAc` does
not match because the already matched Gal still counts for negative
lookbehind.

Lookaround, branch-absence checks, and anchors cannot be repeated.
Assertions alone do not produce a match because every returned match
must contain at least one residue.

### Capturing part of a match

Use `(?<name>...)` to give a matched path a name. For example,
`(?<terminal>Gal)-GlcNAc` matches the two-residue path and captures its
Gal as `terminal`.
[`rex_match()`](https://glycoverse.github.io/glyrex/reference/rex_query.md)
and
[`rex_match_all()`](https://glycoverse.github.io/glyrex/reference/rex_query.md)
return the full match in the `match` column and captured fragments in
their named columns. This named-capture syntax is an extension to
glycowork.

Capture names must start with a letter, contain only letters, digits, or
underscores, and be unique within the pattern. The name `match` is
reserved. Use `(?:...)` when grouping without a named capture is
desired. See
[`rex_detect()`](https://glycoverse.github.io/glyrex/reference/rex_query.md)
for how matches are ordered and counted, including overlapping matches
and matches contained within larger ones.

## Examples

``` r
p <- rex_compile("Galb4-(Fuca3-)GlcNAc")
rex_detect("Gal(b1-4)[Fuc(a1-3)]GlcNAc", p)
#> [1] TRUE
rex_detect(c("GlcNAc(a1-", "GlcNAc(b1-"), "GlcNAca$")
#> [1]  TRUE FALSE
```
