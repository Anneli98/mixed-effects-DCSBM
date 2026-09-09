# ============================================================
# estimate_membership.R
# RJSE community-membership estimation
# ============================================================
#
# One public call performs:
#   1) community-free auxiliary fitting;
#   2) Pearson residual construction;
#   3) residual sum-of-squares aggregation and spectral clustering.
#
# The auxiliary model and spectral workflow are implemented in
# separate internal files under R/subfunctions/community/.
# ============================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. Install it with install.packages('here').",
    call. = FALSE
  )
}

here::i_am("R/estimate_membership.R")

source(
  here::here(
    "R",
    "subfunctions",
    "load_subfunctions.R"
  )
)


#' Estimate shared community memberships by RJSE
#'
#' @param Y An M x E x T edge-only binary network array.
#' @param G Number of communities.
#' @param n Number of nodes. If NULL, inferred from E.
#' @param self_loop Whether self-loops are included.
#' @param n_cores Number of workers for the auxiliary fits. If NULL,
#'   choose automatically.
#' @param seed Optional random seed used by spectral clustering.
#' @param return_auxiliary Whether to return the auxiliary fit for diagnostics.
#' @param verbose Whether to print auxiliary / aggregation progress.
#' @param aux_control Optional named list of auxiliary-model controls.
#'
#' @return An object containing the estimated membership vector and,
#'   optionally, auxiliary-model details.
estimate_membership <- function(
    Y,
    G,
    n = NULL,
    self_loop = FALSE,
    n_cores = NULL,
    seed = NULL,
    return_auxiliary = FALSE,
    verbose = FALSE,
    aux_control = list()
) {

  if (length(dim(Y)) != 3L) {
    stop(
      "Y must be an M x E x T edge-only array.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(G) ||
      length(G) != 1L ||
      !is.finite(G) ||
      G < 2
  ) {
    stop(
      "G must be an integer greater than or equal to 2.",
      call. = FALSE
    )
  }

  G <- as.integer(G)

  default_control <- list(
    max_iter = 300,
    epsilon = 1e-5,
    damp = 0.8,
    warmstart_max_iter = 30,
    warmstart_epsilon = 1e-3,
    edge_chunk_size = 50000L,
    max_step = 5,
    ll_patience = 5L,
    ll_decrease_tol = 1e-10,
    residual_eps = 1e-6
  )

  if (!is.list(aux_control)) {
    stop(
      "aux_control must be a named list.",
      call. = FALSE
    )
  }

  control <- utils::modifyList(
    default_control,
    aux_control
  )

  locked <- c(
    "Y",
    "n",
    "self_loop",
    "n_cores",
    "verbose"
  )

  control[intersect(names(control), locked)] <- NULL

  aux_args <- c(
    list(
      Y = Y,
      n = n,
      self_loop = self_loop,
      n_cores = n_cores,
      verbose = verbose
    ),
    control
  )

  auxiliary <- do.call(
    fit_rjse_auxiliary_model,
    aux_args
  )

  n_used <- auxiliary$settings$n

  spectral <- rjse_spectral_clustering(
    R = auxiliary$R,
    n = n_used,
    G = G,
    self_loop = self_loop,
    seed = seed,
    verbose = verbose
  )

  out <- list(
    membership = spectral$membership,
    settings = list(
      method = "RJSE",
      G = G,
      auxiliary = "community_free_bernoulli",
      residual = "Pearson",
      aggregation = "ss",
      weighted = TRUE,
      row_normalize = TRUE
    )
  )

  if (return_auxiliary) {
    out$auxiliary <- auxiliary
  }

  class(out) <- c(
    "rjse_membership",
    class(out)
  )

  out
}


print.rjse_membership <- function(x, ...) {

  cat("RJSE community-membership estimate\n")
  cat(sprintf("  communities: %d\n", x$settings$G))
  cat(sprintf("  nodes: %d\n", length(x$membership)))

  invisible(x)
}
