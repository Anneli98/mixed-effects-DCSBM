# ============================================================
# rjse_spectral_clustering.R
# Residual aggregation and spectral community recovery
# ============================================================
#
# This file intentionally preserves the networkSoSD workflow used
# in the original RJSE detection script:
#
#   aggregate_networks(..., method = "ss")
#   spectral_clustering(..., weighted = TRUE,
#                       row_normalize = TRUE)
#
# No replacement spectral algorithm is used here.
# ============================================================


# ------------------------------------------------------------
# Convert residual edge arrays into a list of residual matrices
# ------------------------------------------------------------
rjse_residual_matrix_list <- function(
    R,
    n,
    self_loop = FALSE
) {

  if (length(dim(R)) != 3L) {
    stop(
      "R must be an M x E x T residual array.",
      call. = FALSE
    )
  }

  M <- dim(R)[1]
  nT <- dim(R)[3]

  A_residual <- edges_to_adjacency_array(
    R,
    n = n,
    self_loop = self_loop
  )

  out <- vector(
    "list",
    M * nT
  )

  for (t in seq_len(nT)) {
    for (m in seq_len(M)) {
      idx <- m + (t - 1L) * M
      out[[idx]] <- A_residual[m, , , t]
    }
  }

  out
}


# ------------------------------------------------------------
# Original residual sum-of-squares spectral workflow
# ------------------------------------------------------------
rjse_spectral_clustering <- function(
    R,
    n,
    G,
    self_loop = FALSE,
    seed = NULL,
    verbose = FALSE
) {

  if (!requireNamespace("networkSoSD", quietly = TRUE)) {
    stop(
      paste0(
        "Package 'networkSoSD' is required for RJSE community ",
        "estimation."
      ),
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

  A_list_residual <- rjse_residual_matrix_list(
    R = R,
    n = n,
    self_loop = self_loop
  )

  A_agg_residual <- networkSoSD::aggregate_networks(
    A_list_residual,
    method = "ss",
    verbose = verbose
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  membership <- networkSoSD::spectral_clustering(
    A_agg_residual,
    G,
    weighted = TRUE,
    row_normalize = TRUE
  )

  list(
    membership = as.integer(membership),
    aggregate = A_agg_residual
  )
}
