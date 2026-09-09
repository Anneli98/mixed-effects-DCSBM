# ============================================================
# Bandwidth selection for the time-varying mixed-effects DCSBM
# ============================================================
# Common bandwidth selected by leave-one-time-point-out
# predictive marginal log-likelihood.

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. Install it with install.packages('here').",
    call. = FALSE
  )
}

here::i_am("R/bandwidth_selection.R")

source(
  here::here(
    "R",
    "subfunctions",
    "load_subfunctions.R"
  )
)


select_bandwidth <- function(
    Y,
    c,
    h_grid,
    kernel = "epanechnikov",
    n_cores = NULL,
    cv_nAGQ = 10,
    max_iter = 1000,
    Q = 10,
    epsilon = 1e-5,
    damp = 0.8,
    omega_newton_maxit = 20,
    omega_newton_tol = 1e-5,
    keep_best_loo_fit = TRUE,
    verbose = TRUE
) {

  .medcsbm_require("statmod")
  .medcsbm_require("lme4")

  if (length(dim(Y)) != 3L) {
    stop(
      "Y must be an M x E x T edge-only array.",
      call. = FALSE
    )
  }

  h_grid <- sort(
    unique(
      as.numeric(h_grid)
    )
  )

  h_grid <- h_grid[
    is.finite(h_grid) &
      h_grid > 0
  ]

  if (length(h_grid) == 0L) {
    stop(
      "h_grid must contain at least one positive bandwidth.",
      call. = FALSE
    )
  }

  c <- as.integer(
    factor(
      c,
      levels = sort(unique(c))
    )
  )

  n <- length(c)
  M <- dim(Y)[1]
  nT <- dim(Y)[3]
  G <- length(unique(c))

  eg <- prepare_edge_index(
    n,
    self_loop = FALSE
  )

  if (eg$E != dim(Y)[2]) {
    stop(
      "The edge dimension of Y is incompatible with c.",
      call. = FALSE
    )
  }

  ker <- get_kernel(kernel)

  n_workers <- resolve_n_cores(
    n_cores,
    n_tasks = nT
  )

  # Internal numerical starting values, computed at observed times
  # and reused only when the corresponding layer belongs to the
  # LOTO training set.
  .start_values <- function(Y_t) {

    ci <- c[eg$i_idx]
    cj <- c[eg$j_idx]

    a <- pmin(ci, cj)
    b <- pmax(ci, cj)

    block_key <- a + (b - 1L) * G

    trials_tbl <- rowsum(
      rep(1, length(block_key)),
      block_key,
      reorder = FALSE
    )

    key_all <- as.integer(
      rownames(trials_tbl)
    )

    trials_all <- as.numeric(
      trials_tbl
    )

    bb_all <- ((key_all - 1L) %/% G) + 1L
    aa_all <- key_all - (bb_all - 1L) * G

    block_levels <- paste0(
      aa_all,
      "_",
      bb_all
    )

    rows <- vector(
      "list",
      M
    )

    for (m in seq_len(M)) {

      succ_tbl <- rowsum(
        Y_t[m, ],
        block_key,
        reorder = FALSE
      )

      key_s <- as.integer(
        rownames(succ_tbl)
      )

      succ_s <- as.numeric(
        succ_tbl
      )

      success_all <- numeric(
        length(key_all)
      )

      names(success_all) <- key_all

      success_all[
        as.character(key_s)
      ] <- succ_s

      rows[[m]] <- data.frame(
        subject = factor(m),
        block = factor(
          paste0(
            aa_all,
            "_",
            bb_all
          ),
          levels = block_levels
        ),
        success = as.numeric(success_all),
        trial = as.numeric(trials_all)
      )
    }

    dat <- do.call(
      rbind,
      rows
    )

    mod <- lme4::glmer(
      cbind(success, trial - success) ~
        -1 + block + (1 | subject),
      data = dat,
      family = binomial,
      nAGQ = 15,
      control = lme4::glmerControl(
        optCtrl = list(
          maxfun = 5e4
        )
      )
    )

    nm <- names(
      lme4::fixef(mod)
    )

    mm <- regexec(
      "^block(\\d+)_(\\d+)$",
      nm
    )

    rr <- regmatches(
      nm,
      mm
    )

    if (!all(lengths(rr) == 3L)) {
      stop(
        "Unexpected GLMM block coefficient names.",
        call. = FALSE
      )
    }

    ii <- as.integer(
      vapply(
        rr,
        `[[`,
        "",
        2L
      )
    )

    jj <- as.integer(
      vapply(
        rr,
        `[[`,
        "",
        3L
      )
    )

    Z0 <- matrix(
      NA_real_,
      G,
      G
    )

    Z0[cbind(ii, jj)] <-
      as.numeric(
        lme4::fixef(mod)
      )

    Z0[cbind(jj, ii)] <-
      as.numeric(
        lme4::fixef(mod)
      )

    vsub <- as.numeric(
      lme4::VarCorr(mod)$subject
    )

    sigma0 <- if (
      !is.na(vsub) &&
        vsub > 0
    ) {
      sqrt(vsub)
    } else {
      runif(1, 0, 1)
    }

    list(
      theta = numeric(n),
      Z = Z0,
      sigma = sigma0
    )
  }

  if (verbose) {
    message(
      "[CV] precomputing numerical starts"
    )
  }

  init_cache <- lapply(
    seq_len(nT),
    function(tt) {
      .start_values(
        Y[
          ,
          ,
          tt,
          drop = FALSE
        ][, , 1]
      )
    }
  )

  valid_h <- vapply(
    h_grid,
    function(h) {
      tryCatch(
        {
          build_time_neighborhoods(
            nT = nT,
            h = h,
            K = ker$fun,
            smooth = TRUE,
            leave_one_time_out = TRUE
          )
          TRUE
        },
        error = function(e) FALSE
      )
    },
    logical(1)
  )

  results <- data.frame(
    h = h_grid,
    cv = -Inf,
    total_loglik = -Inf,
    status = ifelse(
      valid_h,
      "pending",
      "invalid_loo_neighborhood"
    ),
    stringsAsFactors = FALSE
  )

  best_fit <- NULL
  best_cv <- -Inf
  best_h <- NA_real_

  for (ii in seq_along(h_grid)) {

    h <- h_grid[ii]

    if (!valid_h[ii]) {
      next
    }

    if (verbose) {
      message(
        sprintf(
          "[CV] fitting h=%.6g (%d/%d)",
          h,
          ii,
          length(h_grid)
        )
      )
    }

    tw_list <- build_time_neighborhoods(
      nT = nT,
      h = h,
      K = ker$fun,
      smooth = TRUE,
      leave_one_time_out = TRUE
    )

    theta0 <- matrix(
      NA_real_,
      n,
      nT
    )

    Z0 <- array(
      NA_real_,
      dim = c(G, G, nT)
    )

    sigma0 <- numeric(nT)

    for (tt in seq_len(nT)) {

      t_pick <- tw_list[[tt]]$idx[
        which.min(
          abs(
            tw_list[[tt]]$idx - tt
          )
        )
      ]

      start_t <- init_cache[[t_pick]]

      theta0[, tt] <-
        center_theta_by_community(
          start_t$theta,
          c
        )

      Z0[, , tt] <-
        symmetrize_Z(
          start_t$Z
        )

      sigma0[tt] <-
        max(
          as.numeric(start_t$sigma),
          1e-8
        )
    }

    one <- tryCatch(
      {

        loo_fit <- fit_laem_path(
          Y = Y,
          c = c,
          h = h,
          kernel = ker$fun,
          self_loop = FALSE,
          smooth = TRUE,
          leave_one_time_out = TRUE,
          max_iter = max_iter,
          Q = Q,
          epsilon = epsilon,
          damp = damp,
          verbose = FALSE,

          # Use the same original Q-based convergence behavior.
          compute_Q = TRUE,

          omega_newton_maxit = omega_newton_maxit,
          omega_newton_tol = omega_newton_tol,

          initial_values = list(
            theta = theta0,
            Z = Z0,
            sigma = sigma0
          ),

          n_cores = n_workers
        )

        loo_fit$settings$kernel <- ker$name
        loo_fit$settings$kernel_support <- ker$support

        score <- heldout_marginal_loglik(
          Y,
          c,
          loo_fit,
          nAGQ = cv_nAGQ
        )

        list(
          fit = loo_fit,
          score = score
        )
      },
      error = function(e) e
    )

    if (inherits(one, "error")) {

      results$status[ii] <- paste0(
        "error: ",
        conditionMessage(one)
      )

      if (verbose) {
        message(
          sprintf(
            "[CV] h=%.6g failed: %s",
            h,
            conditionMessage(one)
          )
        )
      }

      next
    }

    results$cv[ii] <- one$score$mean
    results$total_loglik[ii] <-
      one$score$total
    results$status[ii] <- "ok"

    if (
      is.finite(one$score$mean) &&
        one$score$mean > best_cv
    ) {

      best_cv <- one$score$mean
      best_h <- h

      if (keep_best_loo_fit) {
        best_fit <- one$fit
      }
    }

    if (verbose) {
      message(
        sprintf(
          "[CV] h=%.6g, CV=%.6f",
          h,
          one$score$mean
        )
      )
    }
  }

  if (!is.finite(best_cv)) {
    stop(
      "All candidate bandwidths failed during leave-one-time-out CV.",
      call. = FALSE
    )
  }

  out <- list(
    h_opt = best_h,
    cv_opt = best_cv,
    scores = results,
    best_loo_fit = best_fit,
    definition = paste(
      "mean held-out marginal log-likelihood",
      "over M subjects and T time points"
    ),
    M = M,
    T = nT
  )

  class(out) <- c(
    "mixed_effects_bandwidth",
    class(out)
  )

  out
}


print.mixed_effects_bandwidth <- function(x, ...) {

  cat(
    sprintf(
      "Selected bandwidth: h = %.6g\n",
      x$h_opt
    )
  )

  cat(
    sprintf(
      "CV score: %.6f\n",
      x$cv_opt
    )
  )

  invisible(x)
}
