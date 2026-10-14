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

Compilation validates the pattern syntax and prepares it for reuse with
query functions such as
[`rex_detect()`](https://glycoverse.github.io/glyrex/reference/rex_detect.md)
and
[`rex_extract()`](https://glycoverse.github.io/glyrex/reference/rex_extract.md).
Passing an already compiled pattern returns it unchanged.

## See also

[`vignette("glyco-regex-syntax", package = "glyrex")`](https://glycoverse.github.io/glyrex/articles/glyco-regex-syntax.md)
for a step-by-step guide to pattern syntax. See
[`rex_detect()`](https://glycoverse.github.io/glyrex/reference/rex_detect.md)
for matching and result conventions.

## Examples

``` r
p <- rex_compile("Galb4-(Fuca3-)GlcNAc")
rex_detect("Gal(b1-4)[Fuc(a1-3)]GlcNAc", p)
#> [1] TRUE
rex_detect(c("GlcNAc(a1-", "GlcNAc(b1-"), "GlcNAca$")
#> [1]  TRUE FALSE
```
