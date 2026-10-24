# A fixed path needs no backtracking. Build this plan once per query pattern,
# rather than caching it on a user-modifiable compiled object.
rex_prepare_pattern <- function(pattern) {
  pattern$path <- rex_path_plan(pattern$tree)
  if (!is.null(pattern$path)) {
    return(pattern)
  }
  index <- 0L
  prepare <- function(term) {
    if (term$kind == "atom") {
      index <<- index + 1L
      term$atom_id <- as.character(index)
    }
    if (!is.null(term$alternatives)) {
      term$alternatives <- lapply(term$alternatives, function(terms) {
        lapply(terms, prepare)
      })
    }
    if (!is.null(term$child)) {
      term$child <- prepare(term$child)
    }
    term
  }
  pattern$tree <- lapply(pattern$tree, prepare)
  pattern
}

rex_path_plan <- function(tree) {
  atoms <- list()
  captures <- list()
  visit <- function(terms) {
    for (term in terms) {
      if (term$kind == "atom") {
        atoms[[length(atoms) + 1L]] <<- term
      } else if (
        term$kind %in% c("capture", "group") && length(term$alternatives) == 1L
      ) {
        start <- length(atoms) + 1L
        if (!visit(term$alternatives[[1L]])) {
          return(FALSE)
        }
        if (term$kind == "capture") {
          captures[[term$name]] <<- seq.int(start, length(atoms))
        }
      } else {
        return(FALSE)
      }
    }
    TRUE
  }
  if (!visit(tree) || !length(atoms)) {
    return(NULL)
  }
  list(atoms = atoms, captures = captures)
}

rex_atom_mask <- function(term, ctx) {
  hit <- rep(FALSE, ctx$n)
  for (spec in term$specs) {
    label <- spec$mono
    if (label %in% c(".", "Monosaccharide")) {
      hit[] <- TRUE
      break
    }
    hit <- hit |
      ((ctx$mono == label |
        ctx$generic == label |
        (label == "Sia" & ctx$mono %in% c("Neu5Ac", "Neu5Gc", "Kdn"))) &
        ctx$sub == spec$sub)
  }
  if (term$root) {
    hit[-ctx$root] <- FALSE
  }
  if (
    !is.null(term$anomer) || !is.null(term$donor) || !is.null(term$acceptor)
  ) {
    ids <- which(hit)
    hit[ids] <- vapply(ids, function(id) rex_atom(term, id, ctx), logical(1))
  }
  hit
}

rex_search_path <- function(ctx, plan, first = FALSE, groups = NULL) {
  # All paths have the same nonzero length, so no strict-subset pruning is
  # needed. Retain the established minimum-node start key and result order.
  ids <- seq_len(ctx$n)
  paths <- matrix(integer(), nrow = ctx$n, ncol = 0L)
  for (atom in plan$atoms) {
    mask <- c(FALSE, rex_atom_mask(atom, ctx))
    keep <- mask[ids + 1L]
    if (!any(keep)) {
      return(list())
    }
    paths <- cbind(paths[keep, , drop = FALSE], ids[keep])
    ids <- ctx$parent[ids[keep]]
  }
  starts <- paths[, 1L]
  if (ncol(paths) > 1L) {
    for (j in 2:ncol(paths)) {
      starts <- pmin.int(starts, paths[, j])
    }
  }
  keep <- which(!duplicated(starts))
  keep <- keep[order(starts[keep])]
  if (first) {
    keep <- if (is.null(groups)) {
      keep[1L]
    } else {
      keep[!duplicated(groups[paths[keep, 1L]])]
    }
  }
  lapply(keep, function(i) {
    nodes <- unname(paths[i, ])
    list(
      cur = ids[i],
      nodes = nodes,
      captures = lapply(plan$captures, function(at) nodes[at]),
      last = nodes[length(nodes)]
    )
  })
}

# Unconstrained fixed paths can evaluate all graph residues in a single pass.
# Offset parent IDs keep every path inside its original graph; results are
# translated back to the original node IDs before they leave this helper.
rex_search_paths <- function(contexts, plan, first = FALSE) {
  sizes <- vapply(contexts, `[[`, integer(1), "n")
  offsets <- cumsum(sizes) - sizes
  parent <- unlist(lapply(contexts, `[[`, "parent"), use.names = FALSE)
  node_offsets <- rep(offsets, sizes)
  parent[parent != 0L] <- parent[parent != 0L] + node_offsets[parent != 0L]
  ctx <- list(
    n = sum(sizes),
    parent = parent,
    mono = unlist(lapply(contexts, `[[`, "mono"), use.names = FALSE),
    generic = unlist(lapply(contexts, `[[`, "generic"), use.names = FALSE),
    sub = unlist(lapply(contexts, `[[`, "sub"), use.names = FALSE)
  )
  membership <- rep(seq_along(contexts), sizes)
  hits <- rex_search_path(ctx, plan, first = first, groups = membership)
  out <- rep(list(list()), length(contexts))
  for (hit in hits) {
    i <- membership[hit$nodes[1L]]
    offset <- offsets[[i]]
    hit$nodes <- hit$nodes - offset
    hit$last <- hit$last - offset
    if (hit$cur) {
      hit$cur <- hit$cur - offset
    }
    hit$captures <- lapply(hit$captures, function(nodes) nodes - offset)
    out[[i]][[length(out[[i]]) + 1L]] <- hit
  }
  out
}
