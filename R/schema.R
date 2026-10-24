# Only immutable residue metadata is retained between calls. Query results and
# graph contexts are always scoped to a single call.
rex_residues <- local({
  schema <- NULL
  function() {
    if (is.null(schema)) {
      mono <- glyrepr::available_monosaccharides()
      schema <<- list(
        known = c(mono, "Sia", "Monosaccharide", "."),
        generic = structure(glyrepr::convert_to_generic(mono), names = mono)
      )
    }
    schema
  }
})
