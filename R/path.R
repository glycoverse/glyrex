# A fixed path needs no backtracking. Build this plan once per query pattern,
# rather than caching it on a user-modifiable compiled object.
rex_prepare_pattern <- function(pattern) {
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
  pattern$path <- rex_path_plan(pattern$tree)
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

rex_search_path <- function(ctx, plan) {
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
  starts <- apply(paths, 1L, min)
  keep <- which(!duplicated(starts))
  keep <- keep[order(starts[keep])]
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
