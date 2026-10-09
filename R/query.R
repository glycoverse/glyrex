#' Query Glycan Structures with Regular Expressions
#'
#' A stringr-like interface to graph-native glycan pattern matching.
#'
#' @param glycans A character vector of glycan expressions in any format
#'   supported by [glyparse::auto_parse()], a `glyrepr_structure` vector, or
#'   one glyrepr-compatible `igraph`. Character inputs may mix formats;
#'   malformed structures throw an error.
#' @param pattern A character vector of patterns or one [rex_compile()] object.
#'   Inputs recycle only from length one. Missing glycans or patterns propagate.
#' @param negate Invert detection when `TRUE`; missing values remain missing.
#' @returns
#' `rex_detect()` returns logicals; `rex_count()` returns integers.
#' `rex_extract()` returns the first match, or `NA`; `rex_extract_all()`
#' returns a list of vectors (empty for no match, length-one `NA` for missing).
#' Extraction returns canonical glyrepr IUPAC strings for character input and
#' glyrepr structure vectors for structure or graph input. A cut fragment
#' retains the anomer and donor position of its former attachment linkage.
#'
#' `rex_match()` returns a character matrix with `match` followed by named
#' capture columns; unmatched captures are `NA`. `rex_match_all()` returns
#' a list of such matrices, with zero rows for no match and one missing row
#' for missing input. Vectorized patterns must have identical capture names.
#'
#' `rex_locate()` returns a data frame with integer list-columns `nodes` and
#' `edges`; `rex_locate_all()` returns a list of these data frames. IDs are
#' one-based and refer to the input graph, or the canonical graph obtained
#' from character input, as in `glyrepr::structure_nodes()` and
#' `glyrepr::structure_edges()`. Only edges within the match are returned.
#' First locations use `NA_integer_` for no match or missing input; all
#' locations use zero rows for no match and one missing row for missing input.
#'
#' `rex_subset()` preserves the input representation and names;
#' `rex_which()` returns matching indices. Both omit missing results and
#' require `pattern` to have length one or the length of `glycans`.
#'
#' @details
#' The matching unit is a connected structural trace. As in glycowork, the
#' first successful trace per starting node is retained, then strict subsets
#' of other traces are removed. Distinct branches with equal sequences remain
#' distinct matches. Matches are ordered by decreasing residue count and then
#' starting node ID. Repeats are greedy unless marked lazy; successful traces
#' may overlap. Counts therefore need not equal text-regex counts.
#' Only connected rooted trees are supported; floating topology or floating
#' substituents are rejected explicitly. See [rex_compile()] for syntax.
#'
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_detect(x, "Gal-HexNAc")
#' rex_count(x, "HexNAc")
#' rex_extract_all(x, "HexNAc")
#' rex_match(x, "(?<terminal>Gal)-GlcNAc")
#' rex_locate_all(x, "Gal-GlcNAc")
#' rex_subset(x, "Gal")
#' rex_which(x, "Fuc")
#' @name rex_query
NULL

rex_queries <- function(glycans, pattern) {
  graph_input <- inherits(glycans, "igraph")
  structure_input <- inherits(glycans, "glyrepr_structure")
  if (!graph_input && !structure_input && !is.character(glycans)) {
    stop(
      "`glycans` must be character, glyrepr_structure, or a glycan igraph.",
      call. = FALSE
    )
  }
  compiled <- inherits(pattern, "glyrex_pattern")
  if (!compiled && !is.character(pattern)) {
    stop("`pattern` must be character or compiled.", call. = FALSE)
  }
  patterns <- if (compiled) {
    list(pattern)
  } else {
    lapply(pattern, function(p) {
      if (is.na(p)) NULL else rex_compile(p)
    })
  }
  nx <- if (graph_input) 1L else length(glycans)
  np <- length(patterns)
  if (nx != np && nx != 1L && np != 1L) {
    stop("Inputs must have equal lengths or length one.", call. = FALSE)
  }
  n <- if (nx == 0L || np == 0L) 0L else max(nx, np)
  graphs <- if (graph_input) {
    list(glycans)
  } else if (structure_input) {
    as.list(glycans)
  } else {
    as.list(glyparse::auto_parse(glycans, on_failure = "error"))
  }
  schemas <- lapply(patterns, function(p) if (is.null(p)) NULL else p$captures)
  patterns <- rep(patterns, length.out = n)
  indices <- rep(seq_len(nx), length.out = n)
  contexts <- lapply(graphs, function(g) {
    if (inherits(g, "igraph")) rex_graph(g) else NULL
  })
  results <- lapply(seq_len(n), function(i) {
    g <- graphs[[indices[i]]]
    p <- patterns[[i]]
    missing <- is.null(p) || !inherits(g, "igraph")
    ctx <- NULL
    hits <- list()
    if (!missing) {
      ctx <- contexts[[indices[i]]]
      hits <- rex_search(ctx, p)
    }
    list(ctx = ctx, pattern = p, hits = hits, missing = missing)
  })
  if (!graph_input && !is.null(names(glycans))) {
    names(results) <- rep(names(glycans), length.out = n)
  }
  attr(results, "schemas") <- schemas
  attr(results, "structure") <- graph_input || structure_input
  results
}

rex_negate <- function(negate) {
  if (!is.logical(negate) || length(negate) != 1L || is.na(negate)) {
    stop("`negate` must be TRUE or FALSE.", call. = FALSE)
  }
}

#' @rdname rex_query
#' @export
rex_detect <- function(glycans, pattern, negate = FALSE) {
  rex_negate(negate)
  result <- vapply(
    rex_queries(glycans, pattern),
    function(x) {
      if (x$missing) NA else length(x$hits) > 0L
    },
    logical(1)
  )
  if (negate) !result else result
}

#' @rdname rex_query
#' @export
rex_count <- function(glycans, pattern) {
  vapply(
    rex_queries(glycans, pattern),
    function(x) {
      if (x$missing) NA_integer_ else length(x$hits)
    },
    integer(1)
  )
}

rex_extract_results <- function(results, first) {
  output <- lapply(results, function(x) {
    if (x$missing || (first && !length(x$hits))) {
      return(NA_character_)
    }
    hits <- if (first) utils::head(x$hits, 1L) else x$hits
    vapply(hits, function(hit) rex_fragment(x$ctx, hit$nodes), character(1))
  })
  if (first) {
    output <- vapply(output, identity, character(1))
    if (isTRUE(attr(results, "structure"))) {
      output <- glyrepr::as_glycan_structure(output)
    }
  } else if (isTRUE(attr(results, "structure"))) {
    output <- lapply(output, glyrepr::as_glycan_structure)
  }
  output
}

#' @rdname rex_query
#' @export
rex_extract <- function(glycans, pattern) {
  rex_extract_results(rex_queries(glycans, pattern), TRUE)
}

#' @rdname rex_query
#' @export
rex_extract_all <- function(glycans, pattern) {
  rex_extract_results(rex_queries(glycans, pattern), FALSE)
}

rex_match_results <- function(results, first) {
  schemas <- unique(attr(results, "schemas"))
  schemas <- Filter(Negate(is.null), schemas)
  if (length(schemas) > 1L) {
    stop(
      "Vectorized patterns must have identical capture names.",
      call. = FALSE
    )
  }
  columns <- c("match", if (length(schemas)) schemas[[1L]])
  output <- lapply(results, function(x) {
    hits <- if (first) utils::head(x$hits, 1L) else x$hits
    n <- if (x$missing || first) max(1L, length(hits)) else length(hits)
    out <- matrix(
      NA_character_,
      nrow = n,
      ncol = length(columns),
      dimnames = list(NULL, columns)
    )
    for (i in seq_along(hits)) {
      out[i, 1L] <- rex_fragment(x$ctx, hits[[i]]$nodes)
      for (capture in columns[-1L]) {
        out[i, capture] <- rex_fragment(x$ctx, hits[[i]]$captures[[capture]])
      }
    }
    out
  })
  if (first) {
    if (!length(output)) {
      return(matrix(
        character(),
        0L,
        length(columns),
        dimnames = list(NULL, columns)
      ))
    }
    out <- do.call(rbind, unname(output))
    rownames(out) <- names(results)
    out
  } else {
    output
  }
}

#' @rdname rex_query
#' @export
rex_match <- function(glycans, pattern) {
  rex_match_results(rex_queries(glycans, pattern), TRUE)
}

#' @rdname rex_query
#' @export
rex_match_all <- function(glycans, pattern) {
  rex_match_results(rex_queries(glycans, pattern), FALSE)
}

rex_locations <- function(results, first) {
  output <- lapply(results, function(x) {
    hits <- if (first) utils::head(x$hits, 1L) else x$hits
    if (x$missing || (first && !length(hits))) {
      return(data.frame(
        nodes = I(list(NA_integer_)),
        edges = I(list(NA_integer_))
      ))
    }
    nodes <- lapply(hits, function(hit) sort(unique(hit$nodes)))
    edges <- lapply(nodes, function(ids) {
      which(x$ctx$edges[, 1L] %in% ids & x$ctx$edges[, 2L] %in% ids)
    })
    data.frame(nodes = I(nodes), edges = I(edges))
  })
  if (first) {
    if (!length(output)) {
      return(data.frame(nodes = I(list()), edges = I(list())))
    }
    out <- do.call(rbind, unname(output))
    # Data frame row names cannot preserve duplicated or missing input names.
    rownames(out) <- NULL
    out
  } else {
    output
  }
}

#' @rdname rex_query
#' @export
rex_locate <- function(glycans, pattern) {
  rex_locations(rex_queries(glycans, pattern), TRUE)
}

#' @rdname rex_query
#' @export
rex_locate_all <- function(glycans, pattern) {
  rex_locations(rex_queries(glycans, pattern), FALSE)
}

#' @rdname rex_query
#' @export
rex_which <- function(glycans, pattern, negate = FALSE) {
  nx <- if (inherits(glycans, "igraph")) 1L else length(glycans)
  np <- if (inherits(pattern, "glyrex_pattern")) 1L else length(pattern)
  if (!np %in% c(1L, nx)) {
    stop("Subsetting requires one pattern or one per glycan.", call. = FALSE)
  }
  which(rex_detect(glycans, pattern, negate = negate))
}

#' @rdname rex_query
#' @export
rex_subset <- function(glycans, pattern, negate = FALSE) {
  keep <- rex_which(glycans, pattern, negate = negate)
  if (inherits(glycans, "igraph")) {
    if (length(keep)) glycans else NULL
  } else {
    glycans[keep]
  }
}
