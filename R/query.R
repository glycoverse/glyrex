rex_queries <- function(glycans, pattern, first = FALSE) {
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
  # Compile each distinct source once, including vectorized patterns.
  pattern_index <- if (compiled) 1L else match(pattern, unique(pattern))
  patterns <- if (compiled) {
    list(pattern)
  } else {
    lapply(unique(pattern), function(p) {
      if (is.na(p)) NULL else rex_compile(p)
    })
  }
  nx <- if (graph_input) 1L else length(glycans)
  np <- length(pattern_index)
  if (nx != np && nx != 1L && np != 1L) {
    stop("Inputs must have equal lengths or length one.", call. = FALSE)
  }
  n <- if (nx == 0L || np == 0L) 0L else max(nx, np)
  if (graph_input) {
    graphs <- list(glycans)
    graph_index <- 1L
  } else {
    structures <- if (structure_input) {
      glycans
    } else {
      glyparse::auto_parse(glycans, on_failure = "error")
    }
    keys <- as.character(structures)
    distinct <- !duplicated(keys)
    graph_index <- match(keys, keys[distinct])
    # Read unique graphs without copying every repeated structure.
    graphs <- glyrepr::get_structure_graphs(
      structures[distinct],
      return_list = TRUE
    )
  }
  schemas <- lapply(patterns, function(p) if (is.null(p)) NULL else p$captures)
  patterns <- lapply(patterns, function(p) {
    if (is.null(p)) NULL else rex_prepare_pattern(p)
  })
  pi <- rep(pattern_index, length.out = n)
  gi <- rep(graph_index, length.out = n)
  pairs <- paste(gi, pi, sep = ":")
  distinct <- !duplicated(pairs)
  index <- match(pairs, pairs[distinct])
  pi <- pi[distinct]
  gi <- gi[distinct]
  contexts <- lapply(graphs, function(g) {
    if (inherits(g, "igraph")) rex_graph(g, validate = graph_input) else NULL
  })
  matched <- vector("list", length(gi))
  for (j in unique(pi)) {
    p <- patterns[[j]]
    if (is.null(p$path)) {
      next
    }
    batchable <- all(vapply(
      p$path$atoms,
      function(atom) {
        !atom$root &&
          is.null(atom$anomer) &&
          is.null(atom$donor) &&
          is.null(atom$acceptor)
      },
      logical(1)
    ))
    if (!batchable) {
      next
    }
    ids <- which(pi == j & !vapply(contexts[gi], is.null, logical(1)))
    if (length(ids) > 1L) {
      matched[ids] <- rex_search_paths(contexts[gi[ids]], p$path, first = first)
    }
  }
  results <- lapply(seq_along(gi), function(i) {
    ctx <- contexts[[gi[i]]]
    p <- patterns[[pi[i]]]
    missing <- is.null(p) || is.null(ctx)
    hits <- if (missing) {
      list()
    } else if (!is.null(matched[[i]])) {
      matched[[i]]
    } else {
      rex_search(ctx, p, first = first)
    }
    list(ctx = ctx, pattern = p, hits = hits, missing = missing)
  })
  attr(results, "index") <- index
  if (!graph_input && !is.null(names(glycans))) {
    attr(results, "input_names") <- rep(names(glycans), length.out = n)
  }
  attr(results, "schemas") <- schemas
  attr(results, "structure") <- graph_input || structure_input
  results
}

rex_restore <- function(output, results) {
  output <- output[attr(results, "index")]
  names(output) <- attr(results, "input_names")
  output
}

rex_negate <- function(negate) {
  if (!is.logical(negate) || length(negate) != 1L || is.na(negate)) {
    stop("`negate` must be TRUE or FALSE.", call. = FALSE)
  }
}

#' Detect Structural Matches in Glycans
#'
#' Test whether each glycan contains a structural pattern.
#'
#' @param glycans A character vector of glycan expressions in any format
#'   supported by [glyparse::auto_parse()], a `glyrepr_structure` vector, or
#'   one glyrepr-compatible `igraph`. Character inputs may mix formats;
#'   malformed structures throw an error.
#' @param pattern A character vector of patterns or one [rex_compile()] object.
#'   Inputs recycle only from length one. Missing glycans or patterns propagate.
#' @param negate Invert detection when `TRUE`; missing values remain missing.
#' @returns
#' A logical vector. Missing glycans or patterns return `NA`.
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
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_detect(x, "Gal-HexNAc")
#' rex_detect(x, "Gal", negate = TRUE)
#' @export
rex_detect <- function(glycans, pattern, negate = FALSE) {
  rex_negate(negate)
  queries <- rex_queries(glycans, pattern, first = TRUE)
  result <- rex_restore(
    vapply(
      queries,
      function(x) {
        if (x$missing) NA else length(x$hits) > 0L
      },
      logical(1)
    ),
    queries
  )
  if (negate) !result else result
}

#' Count Structural Matches in Glycans
#'
#' Count the structural matches of a pattern in each glycan.
#'
#' @inheritParams rex_detect
#' @returns
#' An integer vector, with zero for no match and `NA_integer_` for
#' missing glycans or patterns.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_count(x, "HexNAc")
#' @export
rex_count <- function(glycans, pattern) {
  queries <- rex_queries(glycans, pattern)
  rex_restore(
    vapply(
      queries,
      function(x) {
        if (x$missing) NA_integer_ else length(x$hits)
      },
      integer(1)
    ),
    queries
  )
}

rex_extract_results <- function(results, first) {
  prepared <- rex_prepare_fragments(
    results,
    first,
    structures = isTRUE(attr(results, "structure"))
  )
  fragments <- prepared$structures
  fragment_keys <- prepared$keys
  output <- lapply(results, function(x) {
    if (x$missing || (first && !length(x$hits))) {
      return(NA_character_)
    }
    hits <- if (first && length(x$hits)) x$hits[1L] else x$hits
    vapply(hits, function(hit) rex_fragment(x$ctx, hit$nodes), character(1))
  })
  if (first) {
    output <- vapply(output, identity, character(1))
  }
  if (!isTRUE(attr(results, "structure"))) {
    return(rex_restore(output, results))
  }
  canonical <- as.character(fragments)
  graphs <- glyrepr::get_structure_graphs(fragments, return_list = TRUE)
  names(graphs) <- canonical
  graphs <- graphs[!duplicated(canonical)]
  make <- function(values) {
    keys <- canonical[match(values, fragment_keys)]
    names(keys) <- names(values)
    used <- unique(keys[!is.na(keys)])
    glyrepr::new_glycan_structure(
      keys,
      if (length(used)) graphs[used] else list()
    )
  }
  if (first) {
    return(make(rex_restore(output, results)))
  }
  # Different input glycans often yield the same fragment vector. Build its
  # structure representation once, without caching across separate queries.
  keys <- vapply(
    output,
    function(values) {
      paste(match(values, fragment_keys), collapse = ",")
    },
    character(1)
  )
  distinct <- !duplicated(keys)
  converted <- lapply(output[distinct], make)
  rex_restore(converted[match(keys, keys[distinct])], results)
}

#' Extract the First Structural Match
#'
#' Extract the first matching fragment from each glycan.
#'
#' @inheritParams rex_detect
#' @returns
#' A character vector of canonical glyrepr IUPAC strings for character
#' input, or a glyrepr structure vector for structure or graph input.
#' No match or missing input returns a missing value. A cut fragment
#' retains the anomer and donor position of its former attachment linkage.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_extract(x, "HexNAc")
#' @export
rex_extract <- function(glycans, pattern) {
  rex_extract_results(rex_queries(glycans, pattern, first = TRUE), TRUE)
}

#' Extract All Structural Matches
#'
#' Extract all matching fragments from each glycan.
#'
#' @inheritParams rex_detect
#' @returns
#' A list of character vectors of canonical glyrepr IUPAC strings for
#' character input, or glyrepr structure vectors for structure or graph
#' input. No match returns an empty vector; missing input returns a
#' length-one missing vector. A cut fragment retains the anomer and donor
#' position of its former attachment linkage.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_extract_all(x, "HexNAc")
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
  rex_prepare_fragments(results, first, captures = TRUE)
  output <- lapply(results, function(x) {
    hits <- if (first && length(x$hits)) x$hits[1L] else x$hits
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
  output <- rex_restore(output, results)
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
    rownames(out) <- attr(results, "input_names")
    out
  } else {
    output
  }
}

#' Extract the First Structural Match and Named Captures
#'
#' Extract the first match and its named capture groups from each glycan.
#'
#' @inheritParams rex_detect
#' @returns
#' A character matrix with one row per input, a `match` column, and
#' columns for named captures. No match or missing input produces a
#' missing row; unmatched captures are `NA`. Fragments are canonical
#' glyrepr IUPAC strings. Vectorized patterns must have identical capture
#' names.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_match(x, "(?<terminal>Gal)-GlcNAc")
#' @export
rex_match <- function(glycans, pattern) {
  rex_match_results(rex_queries(glycans, pattern, first = TRUE), TRUE)
}

#' Extract All Structural Matches and Named Captures
#'
#' Extract all matches and their named capture groups from each glycan.
#'
#' @inheritParams rex_detect
#' @returns
#' A list of character matrices with a `match` column followed by named
#' capture columns. No match produces zero rows; missing input produces
#' one missing row. Unmatched captures are `NA`. Fragments are canonical
#' glyrepr IUPAC strings. Vectorized patterns must have identical capture
#' names.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_match_all(x, "(?<terminal>Gal)-GlcNAc")
#' @export
rex_match_all <- function(glycans, pattern) {
  rex_match_results(rex_queries(glycans, pattern), FALSE)
}

rex_locations <- function(results, first) {
  if (first) {
    nodes <- lapply(results, function(x) {
      if (x$missing || !length(x$hits)) {
        return(NA_integer_)
      }
      sort.int(x$hits[[1L]]$nodes, method = "quick")
    })
    edges <- lapply(seq_along(results), function(i) {
      ids <- nodes[[i]]
      if (anyNA(ids)) {
        return(NA_integer_)
      }
      rex_location_edges(results[[i]]$ctx, ids)
    })
    index <- attr(results, "index")
    return(rex_location_frame(nodes[index], edges[index]))
  }
  output <- lapply(results, function(x) {
    hits <- x$hits
    if (x$missing) {
      return(rex_location_frame(list(NA_integer_), list(NA_integer_)))
    }
    nodes <- lapply(hits, function(hit) sort.int(hit$nodes, method = "quick"))
    edges <- lapply(nodes, function(ids) rex_location_edges(x$ctx, ids))
    rex_location_frame(nodes, edges)
  })
  rex_restore(output, results)
}

# Both columns are already lists of integer vectors of equal length. Avoid
# data.frame's general coercion and recycling for every result glycan.
rex_location_frame <- function(nodes, edges) {
  n <- length(nodes)
  class(nodes) <- "AsIs"
  class(edges) <- "AsIs"
  structure(
    list(nodes = nodes, edges = edges),
    class = "data.frame",
    row.names = if (n) c(NA_integer_, -as.integer(n)) else integer()
  )
}

rex_location_edges <- function(ctx, nodes) {
  if (length(nodes) < 2L) {
    return(integer())
  }
  if (is.null(ctx$cache$edge_ids)) {
    ids <- integer(ctx$n)
    ids[ctx$edges[, 2L]] <- seq_len(nrow(ctx$edges))
    ctx$cache$edge_ids <- ids
  }
  ids <- ctx$cache$edge_ids[nodes[ctx$parent[nodes] %in% nodes]]
  if (is.unsorted(ids)) sort.int(ids, method = "quick") else ids
}

#' Locate the First Structural Match
#'
#' Locate the nodes and edges of the first match in each glycan.
#'
#' @inheritParams rex_detect
#' @returns
#' A data frame with integer list-columns `nodes` and `edges`, and one
#' row per input. No match or missing input uses `NA_integer_` in both
#' columns. IDs are one-based and refer to the input graph, or the
#' canonical graph obtained from character input, as in
#' [glyrepr::structure_nodes()] and [glyrepr::structure_edges()]. Only
#' edges within the match are returned.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_locate(x, "Gal-GlcNAc")
#' @export
rex_locate <- function(glycans, pattern) {
  rex_locations(rex_queries(glycans, pattern, first = TRUE), TRUE)
}

#' Locate All Structural Matches
#'
#' Locate the nodes and edges of all matches in each glycan.
#'
#' @inheritParams rex_detect
#' @returns
#' A list of data frames with integer list-columns `nodes` and `edges`.
#' No match produces zero rows; missing input produces one row with
#' `NA_integer_` in both columns. IDs are one-based and refer to the input
#' graph, or the canonical graph obtained from character input, as in
#' [glyrepr::structure_nodes()] and [glyrepr::structure_edges()]. Only
#' edges within the match are returned.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_locate_all(x, "Gal-GlcNAc")
#' @export
rex_locate_all <- function(glycans, pattern) {
  rex_locations(rex_queries(glycans, pattern), FALSE)
}

#' Find Indices of Glycans with Structural Matches
#'
#' Return the indices of glycans that match a structural pattern.
#'
#' @inheritParams rex_detect
#' @returns
#' An integer vector of matching input indices. Missing results are
#' omitted. The pattern must have length one or the length of `glycans`.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_which(x, "Fuc")
#' rex_which(x, "Gal", negate = TRUE)
#' @export
rex_which <- function(glycans, pattern, negate = FALSE) {
  nx <- if (inherits(glycans, "igraph")) 1L else length(glycans)
  np <- if (inherits(pattern, "glyrex_pattern")) 1L else length(pattern)
  if (!np %in% c(1L, nx)) {
    stop("Subsetting requires one pattern or one per glycan.", call. = FALSE)
  }
  which(rex_detect(glycans, pattern, negate = negate))
}

#' Subset Glycans by Structural Matches
#'
#' Keep glycans that match a structural pattern.
#'
#' @inheritParams rex_detect
#' @returns
#' A subset preserving the input representation and names. Missing
#' results are omitted. For a single graph input, returns the graph on a
#' match and `NULL` otherwise. The pattern must have length one or the
#' length of `glycans`.
#'
#' @inherit rex_detect details
#' @family glycan queries
#' @examples
#' x <- c("Gal(b1-4)GlcNAc", "Fuc(a1-3)GlcNAc", NA)
#' rex_subset(x, "Gal")
#' rex_subset(x, "Gal", negate = TRUE)
#' @export
rex_subset <- function(glycans, pattern, negate = FALSE) {
  keep <- rex_which(glycans, pattern, negate = negate)
  if (inherits(glycans, "igraph")) {
    if (length(keep)) glycans else NULL
  } else {
    glycans[keep]
  }
}
