rex_graph <- function(g, validate = TRUE) {
  # Structure vectors already contain validated graphs. Raw igraph inputs
  # still need full semantic validation before building a query context.
  if (validate) {
    glyrepr::validate_glycan_graph(g)
  }
  vertices <- igraph::vertex_attr(g)
  attributes <- igraph::graph_attr(g)
  n <- length(vertices$mono)
  if (
    !n ||
      (validate && !igraph::is_tree(g, mode = "out")) ||
      length(attributes$floating_parts) ||
      length(attributes$floating_substituents)
  ) {
    stop(
      "Queries require a connected rooted tree without floating metadata.",
      call. = FALSE
    )
  }
  edges <- igraph::as_edgelist(g, names = FALSE)
  storage.mode(edges) <- "integer"
  parent <- integer(n)
  links <- rep("", n)
  if (nrow(edges)) {
    parent[edges[, 2L]] <- edges[, 1L]
    links[edges[, 2L]] <- igraph::edge_attr(g)$linkage
  }
  root <- which(parent == 0L)
  links[root] <- paste0(attributes$anomer, "-")
  mono <- vertices$mono
  list(
    graph = g,
    n = n,
    parent = parent,
    links = links,
    mono = mono,
    generic = unname(rex_residues()$generic[mono]),
    sub = vertices$sub,
    root = root,
    leaf = !seq_len(n) %in% parent,
    edges = edges,
    fragments = new.env(parent = emptyenv()),
    alditol = isTRUE(attributes$alditol),
    cache = new.env(parent = emptyenv())
  )
}

rex_field <- function(want, have) {
  if (is.null(want) || want == "?") {
    return(TRUE)
  }
  # The glycan must be at least as specific as the pattern.
  values <- strsplit(have, "/", fixed = TRUE)[[1L]]
  length(values) > 0L &&
    all(values %in% strsplit(want, "/", fixed = TRUE)[[1L]])
}

rex_atom <- function(term, id, ctx) {
  if (is.null(id) || id == 0L) {
    return(FALSE)
  }
  residue <- any(vapply(
    term$specs,
    function(spec) {
      label <- spec$mono
      hit <- label %in%
        c(".", "Monosaccharide") ||
        label == ctx$mono[id] ||
        label == ctx$generic[id] ||
        (label == "Sia" && ctx$mono[id] %in% c("Neu5Ac", "Neu5Gc", "Kdn"))
      hit && (label %in% c(".", "Monosaccharide") || spec$sub == ctx$sub[id])
    },
    logical(1)
  ))
  if (!residue || (term$root && id != ctx$root)) {
    return(FALSE)
  }
  if (is.null(term$anomer) && is.null(term$donor) && is.null(term$acceptor)) {
    return(TRUE)
  }
  if (is.null(ctx$cache$values)) {
    ctx$cache$values <- regmatches(
      ctx$links,
      regexec("^([ab?])([0-9?]+)-(.*)$", ctx$links)
    )
  }
  fields <- ctx$cache$values[[id]]
  if (!length(fields)) {
    return(FALSE)
  }
  rex_field(term$anomer, fields[2L]) &&
    rex_field(term$donor, fields[3L]) &&
    (is.null(term$acceptor) ||
      (id != ctx$root && rex_field(term$acceptor, fields[4L])))
}

rex_connected <- function(nodes, ctx) {
  length(nodes) > 0L && sum(ctx$parent[nodes] %in% nodes) == length(nodes) - 1L
}

rex_search <- function(ctx, compiled, first = FALSE) {
  if (!is.null(compiled$path)) {
    return(rex_search_path(ctx, compiled$path, first = first))
  }
  atom_masks <- new.env(parent = emptyenv())
  atom <- function(term, id) {
    if (is.null(term$atom_id)) {
      return(rex_atom(term, id, ctx))
    }
    mask <- atom_masks[[term$atom_id]]
    if (is.null(mask)) {
      mask <- rex_atom_mask(term, ctx)
      atom_masks[[term$atom_id]] <- mask
    }
    mask[id]
  }
  initial <- function(cur = NULL, used = integer()) {
    list(cur = cur, nodes = used, captures = list(), last = 0L)
  }
  sequence <- function(terms, states) {
    for (term in terms) {
      next_states <- list()
      for (state in states) {
        next_states <- c(next_states, evaluate(term, state))
      }
      states <- next_states
      if (!length(states)) break
    }
    states
  }
  alternatives <- function(term, state) {
    out <- list()
    for (alt in term$alternatives) {
      out <- c(out, sequence(alt, list(state)))
    }
    out
  }
  attached <- function(term, state) {
    probe <- initial(used = state$nodes)
    hits <- alternatives(term, probe)
    Filter(function(hit) identical(hit$cur, state$cur), hits)
  }
  evaluate <- function(term, state) {
    kind <- term$kind
    if (kind == "atom" || kind == "not") {
      pool <- if (is.null(state$cur)) seq_len(ctx$n) else state$cur
      out <- list()
      for (id in pool) {
        if (!id || id %in% state$nodes) {
          next
        }
        ok <- if (kind == "atom") {
          atom(term, id)
        } else {
          !any(vapply(
            term$alternatives,
            function(alt) atom(alt[[1L]], id),
            logical(1)
          ))
        }
        if (ok) {
          hit <- state
          hit$nodes <- c(state$nodes, id)
          hit$cur <- ctx$parent[id]
          hit$last <- id
          out[[length(out) + 1L]] <- hit
        }
      }
      return(out)
    }
    if (kind == "anchor") {
      if (term$value == "^") {
        if (is.null(state$cur)) {
          return(lapply(which(ctx$leaf), function(id) {
            hit <- state
            hit$cur <- id
            hit
          }))
        }
        ok <- state$cur > 0L && ctx$leaf[state$cur]
      } else if (term$value == "$") {
        ok <- state$last == ctx$root
      } else {
        id <- if (state$last) state$last else state$cur
        if (is.null(id)) {
          return(lapply(
            which(!ctx$leaf & seq_len(ctx$n) != ctx$root),
            function(id) {
              hit <- state
              hit$cur <- id
              hit
            }
          ))
        }
        ok <- id > 0L && id != ctx$root && !ctx$leaf[id]
      }
      return(if (ok) list(state) else list())
    }
    if (kind == "group") {
      return(alternatives(term, state))
    }
    if (kind == "capture") {
      hits <- alternatives(term, state)
      return(lapply(hits, function(hit) {
        hit$captures[[term$name]] <- setdiff(hit$nodes, state$nodes)
        hit
      }))
    }
    if (kind == "branch") {
      if (is.null(state$cur)) {
        return(alternatives(term, state))
      }
      hits <- attached(term, state)
      return(lapply(hits, function(hit) {
        hit$captures <- utils::modifyList(state$captures, hit$captures)
        hit$last <- state$last
        hit
      }))
    }
    if (kind == "absent") {
      if (is.null(state$cur)) {
        out <- list()
        for (id in seq_len(ctx$n)) {
          probe <- state
          probe$cur <- id
          out <- c(out, evaluate(term, probe))
        }
        return(out)
      }
      return(if (!length(attached(term, state))) list(state) else list())
    }
    if (kind == "look") {
      # At an unspecified initial boundary, test every possible starting node.
      if (is.null(state$cur)) {
        out <- list()
        for (id in seq_len(ctx$n)) {
          probe <- state
          probe$cur <- id
          out <- c(out, evaluate(term, probe))
        }
        return(out)
      }
      behind <- startsWith(term$direction, "?<")
      probe <- initial(state$cur)
      hits <- if (behind) attached(term, probe) else alternatives(term, probe)
      positive <- !endsWith(term$direction, "!")
      return(if ((length(hits) > 0L) == positive) list(state) else list())
    }
    if (kind == "repeat") {
      limit <- min(term$max, ctx$n)
      if (term$min > limit) {
        return(list())
      }
      levels <- list(list(state))
      if (limit > 0L) {
        for (i in seq_len(limit)) {
          hits <- list()
          for (prior in levels[[i]]) {
            next_hits <- evaluate(term$child, prior)
            # Never recurse on zero-width results.
            hits <- c(
              hits,
              if (i == 1L) {
                next_hits
              } else {
                Filter(
                  function(hit) length(hit$nodes) > length(prior$nodes),
                  next_hits
                )
              }
            )
          }
          levels[[i + 1L]] <- hits
          if (!length(hits)) break
        }
      }
      counts <- seq.int(term$min, min(limit, length(levels) - 1L))
      counts <- counts[counts >= term$min & counts <= length(levels) - 1L]
      if (!term$lazy) {
        counts <- rev(counts)
      }
      return(unlist(levels[counts + 1L], recursive = FALSE))
    }
    stop("Invalid compiled syntax tree.", call. = FALSE)
  }
  hits <- sequence(compiled$tree, list(initial()))
  hits <- Filter(function(hit) rex_connected(hit$nodes, ctx), hits)
  # As in glycowork, keep the first successful trace for each start node,
  # then remove traces wholly contained in a larger successful trace.
  if (!length(hits)) {
    return(hits)
  }
  starts <- vapply(hits, function(hit) min(hit$nodes), integer(1))
  hits <- hits[!duplicated(starts)]
  keep <- vapply(
    seq_along(hits),
    function(i) {
      !any(vapply(
        seq_along(hits),
        function(j) {
          length(hits[[i]]$nodes) < length(hits[[j]]$nodes) &&
            all(hits[[i]]$nodes %in% hits[[j]]$nodes)
        },
        logical(1)
      ))
    },
    logical(1)
  )
  hits <- hits[keep]
  sizes <- vapply(hits, function(hit) length(hit$nodes), integer(1))
  starts <- vapply(hits, function(hit) min(hit$nodes), integer(1))
  hits <- hits[order(-sizes, starts)]
  if (first) hits[1L] else hits
}

rex_fragment <- function(ctx, nodes) {
  if (!length(nodes)) {
    return(NA_character_)
  }
  if (is.unsorted(nodes)) {
    nodes <- sort.int(nodes, method = "quick")
  }
  key <- paste(nodes, collapse = ",")
  if (exists(key, ctx$fragments, inherits = FALSE)) {
    return(ctx$fragments[[key]])
  }
  root <- nodes[!ctx$parent[nodes] %in% nodes]
  if (length(root) != 1L) {
    stop("A capture must form a connected structure.", call. = FALSE)
  }
  graph <- igraph::induced_subgraph(ctx$graph, nodes)
  graph <- igraph::set_graph_attr(
    graph,
    "anomer",
    sub("-.*$", "", ctx$links[root])
  )
  if (root != ctx$root) {
    graph <- igraph::set_graph_attr(graph, "alditol", FALSE)
  }
  glyrepr::graph_to_iupac(glyrepr::canonicalize_glycan_graph(graph))
}

# Serialize requested subtrees and canonicalize them in one parser batch.
# Branch ordering remains the responsibility of glyrepr. This avoids repeated
# igraph subgraph creation and canonicalization. The
# fragment cache belongs to this query's context and cannot outlive the query.
rex_prepare_fragments <- function(
  results,
  first,
  captures = FALSE,
  structures = FALSE
) {
  records <- list()
  targets <- list()
  preserve <- character()
  canonicalize <- logical()
  for (result in results) {
    ctx <- result$ctx
    hits <- if (first && length(result$hits)) result$hits[1L] else result$hits
    if (length(hits)) {
      if (is.null(ctx$cache$labels)) {
        ctx$cache$labels <- paste0(
          ctx$mono,
          gsub(",", "", ctx$sub, fixed = TRUE)
        )
        ctx$cache$ends <- sub("-.*$", "-", ctx$links)
      }
      ctx$labels <- ctx$cache$labels
      ctx$ends <- ctx$cache$ends
    }
    for (hit in hits) {
      groups <- c(list(hit$nodes), if (captures) hit$captures)
      for (nodes in groups) {
        if (!length(nodes)) {
          next
        }
        if (is.unsorted(nodes)) {
          nodes <- sort.int(nodes, method = "quick")
        }
        key <- paste(nodes, collapse = ",")
        if (exists(key, ctx$fragments, inherits = FALSE)) {
          next
        }
        root <- nodes[!ctx$parent[nodes] %in% nodes]
        if (length(root) != 1L) {
          stop("A capture must form a connected structure.", call. = FALSE)
        }
        children <- nodes[nodes != root]
        linear <- !anyDuplicated(ctx$parent[children])
        if (linear) {
          path <- integer(length(nodes))
          path[1L] <- nodes[!nodes %in% ctx$parent[nodes]]
          if (length(path) > 1L) {
            for (j in 2:length(path)) {
              path[j] <- ctx$parent[path[j - 1L]]
            }
          }
          labels <- ctx$labels[path]
          if (root == ctx$root && ctx$alditol) {
            labels[length(path)] <- paste0(labels[length(path)], "-ol")
          }
          links <- ctx$links[path]
          links[length(path)] <- ctx$ends[root]
          value <- paste0(
            labels,
            "(",
            links,
            c(rep(")", length(path) - 1L), ""),
            collapse = ""
          )
        } else {
          emit <- function(id) {
            children <- nodes[ctx$parent[nodes] == id]
            label <- ctx$labels[id]
            if (!length(children)) {
              return(label)
            }
            arms <- vapply(
              children,
              function(child) {
                paste0(emit(child), "(", ctx$links[child], ")")
              },
              character(1)
            )
            paste0(
              arms[1L],
              paste0(
                if (length(arms) > 1L) {
                  paste0("[", arms[-1L], "]")
                },
                collapse = ""
              ),
              label
            )
          }
          value <- paste0(
            emit(root),
            if (root == ctx$root && ctx$alditol) "-ol",
            "(",
            sub("-.*$", "", ctx$links[root]),
            "-"
          )
        }
        # Equal/ambiguous branch positions can have stable-order ties. A
        # string round trip can reorder these ties, whereas graph extraction
        # must preserve the original serializer's exact output.
        tied <- FALSE
        if (!linear) {
          positions <- sub(".*-", "", ctx$links[children])
          position_groups <- split(positions, ctx$parent[children])
          tied <- any(vapply(
            position_groups,
            function(pos) {
              length(pos) > 1L &&
                (anyDuplicated(pos) > 0L || any(!grepl("^[0-9]+$", pos)))
            },
            logical(1)
          ))
        }
        if (tied) {
          value <- rex_fragment(ctx, nodes)
        }
        records[[length(records) + 1L]] <- value
        preserve[length(records)] <- if (tied) value else NA_character_
        canonicalize[length(records)] <- !linear && !tied
        targets[[length(targets) + 1L]] <- list(
          cache = ctx$fragments,
          key = key
        )
        ctx$fragments[[key]] <- NA_character_
      }
    }
  }
  if (length(records)) {
    fragments <- unlist(records, use.names = FALSE)
    values <- NULL
    if (structures) {
      values <- glyrepr::as_glycan_structure(fragments)
      fragments[canonicalize] <- as.character(values)[canonicalize]
    } else if (any(canonicalize)) {
      fragments[
        canonicalize
      ] <- as.character(glyrepr::as_glycan_structure(fragments[canonicalize]))
    }
    fragments[!is.na(preserve)] <- preserve[!is.na(preserve)]
    for (i in seq_along(targets)) {
      target <- targets[[i]]
      target$cache[[target$key]] <- fragments[[i]]
    }
    return(list(structures = values, keys = fragments))
  }
  list(
    structures = if (structures) {
      glyrepr::as_glycan_structure(character())
    } else {
      NULL
    },
    keys = character()
  )
}
