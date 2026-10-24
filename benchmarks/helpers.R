# Load both implementations into isolated, byte-compiled environments. No
# installed package or input parsing code is modified by these comparisons.
load_glyrex_version <- function(ref = NULL) {
  env <- new.env(parent = baseenv())
  files <- if (is.null(ref)) {
    list.files("R", pattern = "\\.R$", full.names = TRUE)
  } else {
    system2(
      "git",
      c("ls-tree", "-r", "--name-only", shQuote(ref), "R"),
      stdout = TRUE
    )
  }
  for (file in files) {
    if (is.null(ref)) {
      sys.source(file, envir = env)
    } else {
      code <- system2(
        "git",
        c("show", shQuote(paste0(ref, ":", file))),
        stdout = TRUE
      )
      eval(parse(text = code), envir = env)
    }
  }
  for (name in ls(env)) {
    if (is.function(env[[name]])) env[[name]] <- compiler::cmpfun(env[[name]])
  }
  env
}

query_functions <- c(
  "rex_count",
  "rex_detect",
  "rex_extract",
  "rex_extract_all",
  "rex_locate",
  "rex_locate_all",
  "rex_match",
  "rex_match_all",
  "rex_subset",
  "rex_which"
)

normalize_result <- function(value) {
  if (inherits(value, "glyrepr_structure")) {
    return(list(type = "glyrepr_structure", values = as.character(value)))
  }
  if (inherits(value, "igraph")) {
    return(glyrepr::structure_to_iupac(value))
  }
  if (is.list(value) && !is.data.frame(value)) {
    return(lapply(value, normalize_result))
  }
  value
}
