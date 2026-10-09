#' Compile a Glycan Regular Expression
#'
#' Compile a reusable structural query. Patterns run from the non-reducing end
#' towards the reducing end; they are not text regular expressions.
#'
#' @param pattern One non-missing pattern string, or a compiled pattern.
#' @returns A `glyrex_pattern` containing a parsed syntax tree and capture names.
#'
#' @details
#' ## Reading a structural pattern
#'
#' A pattern describes connected residues in a glycan graph.
#' Read a path from left to right, towards the reducing end.
#' Separate adjacent residues with `-`: for example,
#' `Gal-GlcNAc` asks for a Gal attached directly to a GlcNAc. It does not skip
#' intervening residues. Without end anchors, this path can occur within a
#' larger glycan and does not have to describe the whole structure.
#'
#' Start with residue names, then add constraints as needed:
#'
#' * `Gal` matches a galactose residue.
#' * `.` matches any one residue, not any number of residues.
#' * `Hex` and `HexNAc` match members of the corresponding glyrepr classes;
#'   for example, `Hex` matches Gal, Glc, and Man.
#' * `Sia` matches Neu5Ac, Neu5Gc, and Kdn.
#' * `[Gal|Man]` matches either Gal or Man. Alternatives can also be paths,
#'   as in `[Gal-GlcNAc|Gal-GalNAc]`.
#' * `[^Fuc]` matches one residue that is not Fuc. Negated alternatives must
#'   each describe a single residue; `[^Fuc|Xyl]` excludes both Fuc and Xyl.
#'
#' ## Specifying linkages and the reducing end
#'
#' A linkage belongs to the residue on its non-reducing side. For example,
#' `Galb4-GlcNAc` requires Gal to be attached to GlcNAc by a beta1-4 linkage.
#' The equivalent IUPAC form is `Gal(b1-4)-GlcNAc`: `b` is the anomer,
#' `1` the donor position on Gal, and `4` the acceptor position on GlcNAc.
#' Shorthand such as `Galb4` assumes donor position 1; for Neu5Ac, Neu5Gc,
#' Kdn, and `Sia`, shorthand assumes donor position 2.
#'
#' Omitted linkage fields impose no constraint. Thus `Gal-GlcNAc` does not
#' specify a linkage, and `GlcNAca` specifies only the alpha anomer, either
#' at an internal residue or at the reducing end. A specified acceptor
#' position, as in `GlcNAca4`, requires an outgoing glycosidic linkage towards
#' the reducing end and cannot match the reducing-end residue itself.
#'
#' To require reducing-end GlcNAc with an alpha anomer, use `GlcNAca$`.
#' The glyrepr notation `GlcNAc(a1-` additionally specifies anomeric position
#' 1 and explicitly requires the reducing end. Its closed form,
#' `GlcNAc(a1-)`, can also be used inside groups.
#'
#' Unknown linkage fields (`?`) act as wildcards in either the pattern or
#' the input glycan, following glycowork. A specified linkage can therefore
#' match an input whose corresponding field is unknown; a match establishes
#' compatibility, not certainty about that linkage. Ambiguous positions
#' match when their possible values overlap, so acceptor `3/4` is compatible
#' with `4` but not `6`.
#'
#' ## Requiring or excluding branches
#'
#' Ordinary parentheses describe a sibling branch attached to the next
#' residue in the main path. In `Galb4-(Fuca3)-GlcNAc`, both Gal and Fuc
#' attach to the same GlcNAc: Gal by beta1-4 and Fuc by alpha1-3. The Fuc
#' is not inserted between Gal and GlcNAc. Required branch residues are part
#' of the match. Additional branches are allowed unless explicitly excluded.
#' Pattern parentheses serve a different purpose from the square brackets
#' used for branches in an input IUPAC string such as
#' `Gal(b1-4)[Fuc(a1-3)]GlcNAc`.
#'
#' Use `(!...-)` to require the absence of a sibling branch. For example,
#' `Galb4-(!Fuc-)GlcNAc` matches the Gal-GlcNAc path only if the GlcNAc has
#' no remaining Fuc branch with any linkage. `Galb4-(!Fuca3-)GlcNAc` excludes
#' only a branch compatible with alpha1-3 Fuc. The check considers branches
#' not already used by the preceding pattern and can describe several
#' residues, as in `Gal-(!Gal-GlcNAc-)Man`.
#'
#' Keep the attachment dash inside the absence group, immediately before
#' `)`, and write the attachment residue directly after it. Bare `!Fuc`,
#' `[!Fuc]`, and absence groups without the final dash are rejected.
#'
#' Unlike `[^Fuc]`, which matches and includes one non-Fuc residue,
#' `(!Fuc-)` only checks a condition and adds no residues to the match.
#' Consequently, `[^Fuc]-GlcNAc` requires a non-Fuc residue before GlcNAc;
#' `(!Fuc-)GlcNAc` can match GlcNAc alone if it has no Fuc branch.
#'
#' ## Grouping and repetition
#'
#' Use square brackets or `(?:...)` to group a path. Ordinary parentheses
#' mean a branch, not a path group. Attach a repetition operator to the group:
#'
#' * `[Gal]?`: zero or one Gal.
#' * `[Gal]*`: zero or more consecutive Gal residues.
#' * `[Gal]+`: one or more consecutive Gal residues.
#' * `[Gal]{2}`: exactly two consecutive Gal residues.
#' * `[Gal]{2,4}`: between two and four consecutive Gal residues.
#' * `[Gal]{2,}`: at least two consecutive Gal residues.
#' * `[Gal]{,4}`: zero to four consecutive Gal residues.
#'
#' A group can contain more than one residue: `[Galb4-GlcNAcb3]{2}` repeats
#' the two-residue path twice. Use grouped syntax for repetition rather than
#' forms such as `Gal*`, which are rejected.
#'
#' Repetition is greedy by default: longer repetitions are tried first.
#' Append `?` to try shorter repetitions first, as in `[Gal]+?`. In either
#' case, the rest of the pattern must also match. Open upper bounds are
#' limited by the input graph, not a fixed repeat count. A pattern that
#' permits zero repetitions never produces an empty match.
#'
#' ## Anchors and lookaround assertions
#'
#' Assertions restrict where a match is allowed without adding residues:
#'
#' * `^Gal` requires Gal at a non-reducing terminus (a residue with no
#'   further residues attached on its non-reducing side).
#' * `GlcNAc$` requires GlcNAc at the reducing end.
#' * `Gal%` requires Gal to be internal: neither a non-reducing terminus
#'   nor the reducing-end residue.
#'
#' Lookaround checks a neighbouring path without including it in the result.
#' Lookahead examines the reducing-end direction: `Gal-(?=GlcNAc)` matches
#' Gal only when the next residue is GlcNAc, whereas `Gal-(?!GlcNAc)` requires
#' that the next residue does not match GlcNAc (or is absent). Lookbehind
#' examines the non-reducing direction: `(?<=Gal-)GlcNAc` matches GlcNAc
#' when a Gal is attached to it, and `(?<!Gal-)GlcNAc` requires that no such
#' Gal path exists. Paths inside lookbehind are still written towards the
#' reducing end, just like the rest of the pattern.
#'
#' Lookaround, branch-absence checks, and anchors cannot be repeated.
#' Assertions alone do not produce a match because every returned match
#' must contain at least one residue.
#'
#' ## Capturing part of a match
#'
#' Use `(?<name>...)` to give a matched path a name. For example,
#' `(?<terminal>Gal)-GlcNAc` matches the two-residue path and captures its Gal
#' as `terminal`. [rex_match()] and [rex_match_all()] return the full match
#' in the `match` column and captured fragments in their named columns.
#' This named-capture syntax is an extension to glycowork.
#'
#' Capture names must start with a letter, contain only letters, digits,
#' or underscores, and be unique within the pattern. The name `match` is
#' reserved. Use `(?:...)` when grouping without a named capture is desired.
#' See [rex_detect()] for how matches are ordered and counted, including
#' overlapping matches and matches contained within larger ones.
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
        neg <- peek() == "^"
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
          term <- node("not", alternatives = alternatives)
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
        alternatives <- parse_alts(")", branch_absence = kind == "absent")
        term <- node(
          kind,
          alternatives = alternatives,
          name = name,
          direction = prefix
        )
      } else {
        if (peek() == "!") {
          fail(
            "Bare !residue is not supported; use (!branch-) for branch absence or [^residue] for residue negation."
          )
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
  parse_alts <- function(close, branch_absence = FALSE) {
    parse_alternative <- function() {
      if (close == "]" && peek() == "!") {
        fail(
          "'[!...]' is not supported; use '[^...]' for residue negation or '(!...-)' for branch absence."
        )
      }
      alternative <- parse_seq(close)
      if (branch_absence && substr(pattern, pos - 1L, pos - 1L) != "-") {
        fail(
          "Branch absence requires a trailing attachment dash: use (!branch-)residue."
        )
      }
      alternative
    }
    alternatives <- list(parse_alternative())
    while (peek() == "|") {
      pos <<- pos + 1L
      alternatives[[length(alternatives) + 1L]] <- parse_alternative()
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
