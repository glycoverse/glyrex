#' Compile a Glycan Regular Expression
#'
#' Compile a reusable structural query. Patterns run from the non-reducing end
#' towards the reducing end; they are not text regular expressions.
#'
#' @param pattern One non-missing pattern string, or a compiled pattern.
#' @returns A `glyrex_pattern` containing a parsed syntax tree and capture names.
#' @details
#' Residues are separated by `-`. Linkages accept shorthand (`Galb4`) or
#' IUPAC notation (`Gal(b1-4)`). `.` matches any residue; `Hex`, `HexNAc`,
#' and other glyrepr classes match their concrete members. `Sia` matches
#' Neu5Ac, Neu5Gc, and Kdn. `[Gal|Man]` expresses alternatives; `[^Fuc]`
#' (or `[!Fuc]`) excludes residues. Parentheses describe a sibling branch,
#' as in `Galb4-(Fuca3)-GlcNAc`. A negated linked residue, `!Fuca3`,
#' asserts that such a branch is absent at the current attachment point.
#' `(!...)` always asserts branch absence: `(!Fuc)` forbids a Fuc branch
#' with any linkage, whereas `(!Fuca3)` forbids an alpha1-3 Fuc branch.
#' The assertion checks unconsumed branches attached to the next residue,
#' consumes no nodes, and can contain a multi-residue branch pattern.
#'
#' Bracketed groups accept `?`, `*`, `+`, `{n}`, `{n,m}`, `{n,}`, or
#' `{,m}`. Append `?` for lazy repetition. Open bounds are limited by
#' the input graph, not a fixed repeat count. `^` requires a non-reducing
#' terminus, `$` the reducing end, and `%` an internal residue.
#' Lookahead `(?=...)` / `(?!...)` tests the reducing-end direction;
#' lookbehind `(?<=...)` / `(?<!...)` tests the non-reducing direction.
#' Assertions consume no residues. Empty matches are not returned.
#'
#' `(?<name>...)` captures a path (an extension to glycowork); `(?:...)`
#' groups a path without capture. Ordinary parentheses remain branches.
#' The capture name `match` is reserved for the full-match result column.
#'
#' An omitted anomer is unconstrained. `GlcNAca` constrains the anomer
#' at either an internal residue or the reducing end. `GlcNAc(a1-`
#' explicitly requires the reducing end and its anomeric position, using
#' glyrepr notation. A closed form, `GlcNAc(a1-)`, is also accepted inside
#' groups. A specified acceptor position, e.g. `GlcNAca4`, requires an
#' outgoing glycosidic linkage and cannot match the reducing end.
#' Unknown linkage fields are wildcards on either side, following glycowork;
#' ambiguous positions match when their possible values overlap.
#'
#' @examples
#' p <- rex_compile("Galb4-(Fuca3)-GlcNAc")
#' rex_detect("Gal(b1-4)[Fuc(a1-3)]GlcNAc", p)
#' rex_detect(c("GlcNAc(a1-", "GlcNAc(b1-"), "GlcNAca$")
#' @export
rex_compile <- function(pattern) {
  if (inherits(pattern, "glyrex_pattern")) {
    return(pattern)
  }
  if (
    !is.character(pattern) ||
      length(pattern) != 1L ||
      is.na(pattern) ||
      !nzchar(trimws(pattern))
  ) {
    stop("`pattern` must be one non-empty, non-missing string.", call. = FALSE)
  }
  source <- pattern
  pattern <- gsub("\\s+", "", pattern, perl = TRUE)
  if (startsWith(pattern, "r")) {
    pattern <- substring(pattern, 2L)
  }
  pos <- 1L
  size <- nchar(pattern)
  captures <- character()
  peek <- function() if (pos <= size) substr(pattern, pos, pos) else ""
  tail <- function() substring(pattern, pos)
  fail <- function(message) {
    stop(
      sprintf("Invalid pattern at character %d: %s", pos, message),
      call. = FALSE
    )
  }
  take <- function(rx) {
    m <- regexpr(rx, tail(), perl = TRUE)
    if (m[1L] != 1L) {
      return(NULL)
    }
    value <- regmatches(tail(), m)
    pos <<- pos + nchar(value)
    value
  }
  node <- function(kind, ...) list(kind = kind, ...)
  parse_seq <- function(close = "") {
    out <- list()
    while (pos <= size && !peek() %in% c(close, "|")) {
      if (peek() == "-") {
        pos <<- pos + 1L
        if ((!length(out) && close == "") || peek() %in% c("", "|", "-")) {
          fail("Misplaced '-'.")
        }
        next
      }
      ch <- peek()
      grouped <- FALSE
      if (ch %in% c("^", "$", "%")) {
        pos <<- pos + 1L
        term <- node("anchor", value = ch)
      } else if (ch == "[") {
        grouped <- TRUE
        pos <<- pos + 1L
        neg <- peek() %in% c("^", "!")
        if (neg) {
          pos <<- pos + 1L
        }
        alternatives <- parse_alts("]")
        if (neg) {
          if (
            any(vapply(
              alternatives,
              function(x) length(x) != 1L || x[[1]]$kind != "atom",
              logical(1)
            ))
          ) {
            fail("Negated alternatives must be single residues.")
          }
          term <- node(
            if (
              length(alternatives) == 1L &&
                !is.null(alternatives[[1L]][[1L]]$acceptor)
            ) {
              "absent"
            } else {
              "not"
            },
            alternatives = alternatives
          )
        } else {
          term <- node("group", alternatives = alternatives)
        }
      } else if (ch == "(") {
        grouped <- TRUE
        pos <<- pos + 1L
        kind <- "branch"
        name <- NULL
        prefix <- take("^\\?(?:<=|<!|=|!)")
        if (peek() == "!" && is.null(prefix)) {
          pos <<- pos + 1L
          kind <- "absent"
        } else if (!is.null(prefix)) {
          kind <- "look"
        } else if (!is.null(take("^\\?:"))) {
          kind <- "group"
        } else if (!is.null(take("^\\?<"))) {
          name <- take("^[A-Za-z][A-Za-z0-9_]*")
          if (
            is.null(name) || peek() != ">" || name %in% c(captures, "match")
          ) {
            fail("Invalid, duplicate, or reserved capture name.")
          }
          pos <<- pos + 1L
          captures <<- c(captures, name)
          kind <- "capture"
        } else if (peek() == "?") {
          fail("Unknown group syntax.")
        }
        alternatives <- parse_alts(")")
        term <- node(
          kind,
          alternatives = alternatives,
          name = name,
          direction = prefix
        )
      } else {
        absent <- peek() == "!"
        if (absent) {
          pos <<- pos + 1L
        }
        token <- take("^(?:[DL]-)?[A-Za-z0-9.][A-Za-z0-9./]*")
        if (is.null(token)) {
          fail("Expected a residue or group.")
        }
        # Split shorthand only when the whole token is not a known residue.
        anomer <- donor <- acceptor <- NULL
        root <- FALSE
        known <- c(
          glyrepr::available_monosaccharides(),
          "Sia",
          "Monosaccharide",
          "."
        )
        if (!token %in% known) {
          m <- regexec(
            "^(.*?)([ab])([0-9]+(?:/[0-9]+)*|\\?)?$",
            token,
            perl = TRUE
          )
          parts <- regmatches(token, m)[[1L]]
          if (length(parts)) {
            token <- parts[2L]
            anomer <- parts[3L]
            if (nzchar(parts[4L])) acceptor <- parts[4L]
          }
        }
        if (peek() == "?" && is.null(anomer)) {
          suffix <- take("^\\?(?:[0-9]+(?:/[0-9]+)*|\\?)")
          if (!is.null(suffix)) {
            anomer <- "?"
            acceptor <- substring(suffix, 2L)
          }
        } else if (!is.null(anomer) && is.null(acceptor) && peek() == "?") {
          pos <<- pos + 1L
          acceptor <- "?"
        }
        link <- take("^\\([ab?](?:[0-9]+|\\?)-(?:[0-9]+(?:/[0-9]+)*|\\?)?\\)")
        if (is.null(link)) {
          link <- take("^\\([ab?](?:[0-9]+|\\?)-$")
        }
        if (!is.null(link)) {
          fields <- regmatches(
            link,
            regexec("^\\(([ab?])([0-9]+|\\?)-([^)]*)", link, perl = TRUE)
          )[[1L]]
          anomer <- fields[2L]
          donor <- fields[3L]
          acceptor <- fields[4L]
          root <- !nzchar(acceptor)
          if (root) acceptor <- NULL
        } else if (!is.null(acceptor)) {
          donor <- if (token %in% c("Neu5Ac", "Neu5Gc", "Kdn", "Sia")) {
            "2"
          } else {
            "1"
          }
        }
        residues <- strsplit(token, "/", fixed = TRUE)[[1L]]
        # Let glyrepr parse substituted residue names once, during compilation.
        specs <- lapply(residues, function(residue) {
          if (residue %in% known) {
            return(list(mono = residue, sub = ""))
          }
          parsed <- tryCatch(
            glyrepr::as_glycan_structure(residue),
            error = function(e) NULL
          )
          if (is.null(parsed)) {
            fail(paste0("Unknown residue '", residue, "'."))
          }
          g <- as.list(parsed)[[1L]]
          if (igraph::vcount(g) != 1L) {
            fail("Expected one residue.")
          }
          list(
            mono = igraph::vertex_attr(g, "mono"),
            sub = igraph::vertex_attr(g, "sub")
          )
        })
        term <- node(
          "atom",
          specs = specs,
          anomer = anomer,
          donor = donor,
          acceptor = acceptor,
          root = root
        )
        if (absent) {
          term <- if (!is.null(acceptor) || peek() == "-") {
            node("absent", alternatives = list(list(term)))
          } else {
            node("not", alternatives = list(list(term)))
          }
        }
      }
      quant <- take("^(?:\\{[^}]*\\}|[?*+])")
      if (!is.null(quant)) {
        if (!grouped && quant != "?") {
          fail("Attach repetition to a bracketed group.")
        }
        if (term$kind %in% c("anchor", "look", "absent")) {
          fail("Assertions cannot be repeated.")
        }
        bounds <- switch(
          quant,
          "?" = c(0, 1),
          "*" = c(0, Inf),
          "+" = c(1, Inf),
          NULL
        )
        if (is.null(bounds)) {
          body <- substr(quant, 2L, nchar(quant) - 1L)
          if (
            !grepl("^(?:[0-9]+|[0-9]*,[0-9]*)$", body, perl = TRUE) ||
              body == ","
          ) {
            fail("Invalid repetition range.")
          }
          if (!grepl(",", body, fixed = TRUE)) {
            bounds <- rep(as.numeric(body), 2L)
          } else {
            parts <- strsplit(paste0(body, "x"), ",", fixed = TRUE)[[1L]]
            parts[2L] <- sub("x$", "", parts[2L])
            bounds <- c(
              if (nzchar(parts[1L])) as.numeric(parts[1L]) else 0,
              if (nzchar(parts[2L])) as.numeric(parts[2L]) else Inf
            )
          }
          if (!is.finite(bounds[1L]) || bounds[1L] > bounds[2L]) {
            fail("Invalid repetition bounds.")
          }
        }
        lazy <- peek() == "?"
        if (lazy) {
          pos <<- pos + 1L
        }
        term <- node(
          "repeat",
          child = term,
          min = bounds[1L],
          max = bounds[2L],
          lazy = lazy
        )
      }
      # A lookbehind immediately before a branch constrains that branch,
      # not the shared attachment residue.
      branch <- if (term$kind == "repeat") term$child else term
      if (branch$kind == "branch" && length(out)) {
        pending <- list()
        while (length(out) && out[[length(out)]]$kind == "look") {
          pending <- c(list(out[[length(out)]]), pending)
          out <- utils::head(out, -1L)
        }
        if (length(pending)) {
          branch$alternatives <- lapply(branch$alternatives, function(alt) {
            c(pending, alt)
          })
          if (term$kind == "repeat") term$child <- branch else term <- branch
        }
      }
      out[[length(out) + 1L]] <- term
    }
    if (!length(out)) {
      fail("Empty pattern or alternative.")
    }
    out
  }
  parse_alts <- function(close) {
    alternatives <- list(parse_seq(close))
    while (peek() == "|") {
      pos <<- pos + 1L
      alternatives[[length(alternatives) + 1L]] <- parse_seq(close)
    }
    if (peek() != close) {
      fail(paste0("Expected '", close, "'."))
    }
    pos <<- pos + 1L
    alternatives
  }
  tree <- parse_seq()
  if (pos <= size) {
    fail("Unexpected closing delimiter or ungrouped alternative.")
  }
  structure(
    list(pattern = source, tree = tree, captures = captures),
    class = "glyrex_pattern"
  )
}

#' @export
print.glyrex_pattern <- function(x, ...) {
  cat("<glyrex_pattern> ", x$pattern, "\n", sep = "")
  if (length(x$captures)) {
    cat("Captures: ", paste(x$captures, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
