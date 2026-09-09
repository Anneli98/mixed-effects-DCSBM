# ============================================================
# laem_core.R
# LA-EM point-estimation module
# ============================================================
#
# Contains:
#   - LA-E step
#   - constrained LA-M updates
#   - one-target-time fit_laem_time()
#   - full trajectory fit_laem_path()
#
# Convergence follows the original implementation:
#   1) parameter stability;
#   2) early stop after five consecutive decreases of Q;
#   3) return the parameter values associated with the best recorded Q.
#
# Depends on: base_functions.R
# ============================================================


# ------------------------------------------------------------
# Aggregate edge-wise score contributions
# ------------------------------------------------------------
score_from_residuals <- function(R, i_idx, j_idx, n, ab_idx, G) {

  g_theta <- numeric(n)

  tmp_i <- rowsum(R, i_idx, reorder = FALSE)
  tmp_j <- rowsum(R, j_idx, reorder = FALSE)

  ii <- as.integer(rownames(tmp_i))
  jj <- as.integer(rownames(tmp_j))

  g_theta[ii] <- g_theta[ii] + as.numeric(tmp_i)
  g_theta[jj] <- g_theta[jj] + as.numeric(tmp_j)

  gZ_vec <- numeric(G * G)
  agg <- rowsum(R, ab_idx, reorder = FALSE)
  gZ_vec[as.integer(rownames(agg))] <- as.numeric(agg)

  g_Z <- matrix(
    gZ_vec,
    nrow = G,
    ncol = G,
    byrow = FALSE
  )

  g_Z <- symmetrize_Z(g_Z)

  list(
    g_theta = g_theta,
    g_Z = g_Z
  )
}


# ------------------------------------------------------------
# Diagonal information used by the one-step M-updates
# ------------------------------------------------------------
diagonal_information_from_weights <- function(
    W,
    i_idx,
    j_idx,
    n,
    ab_idx,
    G
) {

  H_theta <- numeric(n)

  tmp_i <- rowsum(W, i_idx, reorder = FALSE)
  tmp_j <- rowsum(W, j_idx, reorder = FALSE)

  ii <- as.integer(rownames(tmp_i))
  jj <- as.integer(rownames(tmp_j))

  H_theta[ii] <- H_theta[ii] + as.numeric(tmp_i)
  H_theta[jj] <- H_theta[jj] + as.numeric(tmp_j)

  HZ_vec <- numeric(G * G)
  agg <- rowsum(W, ab_idx, reorder = FALSE)
  HZ_vec[as.integer(rownames(agg))] <- as.numeric(agg)

  H_Z <- matrix(
    HZ_vec,
    nrow = G,
    ncol = G,
    byrow = FALSE
  )

  H_Z <- symmetrize_Z(H_Z)

  list(
    H_theta = H_theta,
    H_Z = H_Z
  )
}


# ------------------------------------------------------------
# KKT-constrained theta update
# ------------------------------------------------------------
theta_constrained_step <- function(
    S_theta,
    H_pos_theta,
    c,
    damp = 1.0,
    eps = 1e-15
) {

  H_pos_theta <- pmax(H_pos_theta, eps)
  invH <- 1 / H_pos_theta
  b <- damp * invH * S_theta

  grp <- as.integer(
    factor(c, levels = sort(unique(c)))
  )

  G <- max(grp)

  rs_b <- rowsum(b, grp, reorder = FALSE)
  rs_h <- rowsum(invH, grp, reorder = FALSE)

  sum_b <- numeric(G)
  sum_h <- numeric(G)

  sum_b[as.integer(rownames(rs_b))] <- as.numeric(rs_b)
  sum_h[as.integer(rownames(rs_h))] <- as.numeric(rs_h)

  alpha <- sum_b / pmax(sum_h, eps)

  b - invH * alpha[grp]
}


# ------------------------------------------------------------
# LA-E step
# ------------------------------------------------------------
laem_e_step <- function(
    theta,
    Z,
    sigma,
    y_list,
    w_time,
    base_eta,
    gh,
    compute_Q = TRUE,
    omega_newton_maxit = 20,
    omega_newton_tol = 1e-6
) {

  M <- nrow(y_list[[1]])
  E <- length(base_eta)

  x <- gh$nodes
  wts <- gh$weights
  Qn <- length(x)

  R <- numeric(E)
  W <- numeric(E)
  Qval <- 0
  Ew2_bar <- 0

  mu_bar <- numeric(E)
  w_bar <- numeric(E)

  for (s in seq_along(y_list)) {

    y_mat_s <- y_list[[s]]
    ws <- w_time[s]

    for (m in seq_len(M)) {

      yk <- y_mat_s[m, ]

      pm <- posterior_mode_variance(
        y_k = yk,
        base_eta = base_eta,
        sigma = sigma,
        maxit = omega_newton_maxit,
        tol = omega_newton_tol
      )

      mk <- pm$mean
      vk <- pm$var

      Ew2_bar <- Ew2_bar + ws * (mk^2 + vk) / M

      mu_bar[] <- 0
      w_bar[] <- 0

      sv <- sqrt(max(vk, 1e-12))

      if (compute_Q) {

        elog1p <- 0

        for (q in seq_len(Qn)) {

          etaq <- base_eta + (mk + sv * x[q])
          muq <- logistic_prob(etaq)

          mu_bar <- mu_bar + wts[q] * muq
          w_bar <- w_bar + wts[q] * (muq * (1 - muq))

          elog1p <- elog1p +
            wts[q] * sum(log1pexp_stable(etaq))
        }

        Qval <- Qval +
          ws * (sum(yk * (base_eta + mk)) - elog1p) -
          ws * 0.5 * (
            log(sigma^2) +
              (mk^2 + vk) / sigma^2
          )

      } else {

        for (q in seq_len(Qn)) {

          etaq <- base_eta + (mk + sv * x[q])
          muq <- logistic_prob(etaq)

          mu_bar <- mu_bar + wts[q] * muq
          w_bar <- w_bar + wts[q] * (muq * (1 - muq))
        }
      }

      R <- R + ws * (yk - mu_bar)
      W <- W + ws * w_bar
    }
  }

  list(
    R = R,
    W = W,
    Ew2_bar = Ew2_bar,
    Q = Qval
  )
}


# ============================================================
# LA-EM fit at one target time
# ============================================================
fit_laem_time <- function(
    Y_sub,
    c,
    t_global,
    tw,
    eg,
    theta = NULL,
    Z = NULL,
    sigma = NULL,
    max_iter = 1000,
    Q = 10,
    epsilon = 1e-5,
    damp = 0.8,
    verbose = FALSE,
    compute_Q = TRUE,
    omega_newton_maxit = 20,
    omega_newton_tol = 1e-5,
    q_patience = 5L,
    q_decrease_tol = 1e-10
) {

  .medcsbm_require("statmod")

  stopifnot(length(dim(Y_sub)) == 3L)

  E <- dim(Y_sub)[2]
  L <- dim(Y_sub)[3]

  stopifnot(
    length(tw$idx) == L,
    length(tw$w) == L
  )

  n <- length(c)
  c <- as.integer(c)
  G <- length(unique(c))

  i_idx <- eg$i_idx
  j_idx <- eg$j_idx

  stopifnot(
    length(i_idx) == E,
    length(j_idx) == E
  )

  y_list <- lapply(
    seq_len(L),
    function(s) {
      Y_sub[, , s, drop = FALSE][, , 1]
    }
  )

  w_s <- tw$w

  if (is.null(theta)) theta <- numeric(n)
  if (is.null(Z)) Z <- matrix(0, G, G)
  if (is.null(sigma)) sigma <- 0.5

  # Put the starting value in the feasible space once.
  # Subsequent theta updates are constrained through KKT only.
  theta <- center_theta_by_community(theta, c)
  Z <- symmetrize_Z(Z)
  sigma <- max(as.numeric(sigma), 1e-8)

  ab_idx <- make_block_index(
    c,
    i_idx,
    j_idx,
    G
  )

  gh <- statmod::gauss.quad.prob(
    Q,
    "normal"
  )

  Q_path <- if (compute_Q) {
    rep(NA_real_, max_iter)
  } else {
    NULL
  }

  best_Q <- -Inf
  best_iter <- NA_integer_

  best_theta <- theta
  best_Z <- Z
  best_sigma <- sigma

  dec_count <- 0L
  stop_reason <- "max_iter"

  for (iter in seq_len(max_iter)) {

    base_eta <- compute_edge_linear_predictor(
      theta,
      Z,
      c,
      i_idx,
      j_idx
    )

    Est <- laem_e_step(
      theta = theta,
      Z = Z,
      sigma = sigma,
      y_list = y_list,
      w_time = w_s,
      base_eta = base_eta,
      gh = gh,
      compute_Q = compute_Q,
      omega_newton_maxit = omega_newton_maxit,
      omega_newton_tol = omega_newton_tol
    )

    R <- Est$R
    W <- Est$W

    # sigma^2 closed-form update
    sigma2_new <- max(
      Est$Ew2_bar,
      1e-15
    )

    sigma_new <- sqrt(sigma2_new)

    gr <- score_from_residuals(
      R,
      i_idx,
      j_idx,
      n,
      ab_idx,
      G
    )

    Hz <- diagonal_information_from_weights(
      W,
      i_idx,
      j_idx,
      n,
      ab_idx,
      G
    )

    # theta: KKT-constrained one-step update
    delta_theta <- theta_constrained_step(
      S_theta = gr$g_theta,
      H_pos_theta = Hz$H_theta,
      c = c,
      damp = damp
    )

    theta_new <- theta + delta_theta

    # Z: diagonal one-step update + symmetry
    step_Z <- damp * (
      gr$g_Z /
        pmax(Hz$H_Z, 1e-15)
    )

    Z_new <- symmetrize_Z(
      Z + step_Z
    )

    relchg <- max(
      abs(
        c(
          theta_new - theta,
          as.vector(Z_new - Z),
          sigma_new - sigma
        )
      )
    ) /
      (
        1 +
          max(
            abs(
              c(
                theta,
                as.vector(Z),
                sigma
              )
            )
          )
      )

    # Update the current state first, matching the original code.
    theta <- theta_new
    Z <- Z_new
    sigma <- sigma_new

    # Q bookkeeping and best-Q tracking:
    # this ordering intentionally follows the original implementation.
    if (compute_Q) {

      Q_path[iter] <- Est$Q

      if (
        is.finite(Q_path[iter]) &&
          Q_path[iter] > best_Q
      ) {

        best_Q <- Q_path[iter]
        best_iter <- iter

        best_theta <- theta
        best_Z <- Z
        best_sigma <- sigma
      }

      if (iter == 1L) {

        dec_count <- 0L

      } else if (
        Q_path[iter] <
          Q_path[iter - 1L] -
          q_decrease_tol
      ) {

        dec_count <- dec_count + 1L

      } else {

        dec_count <- 0L
      }
    }

    if (verbose) {

      if (compute_Q) {

        cat(
          sprintf(
            paste0(
              "t=%d iter=%d Q=%.6f ",
              "bestQ=%.6f(best@%s) ",
              "sigma=%.6f relchg=%.3e ",
              "dec=%d |S_t|=%d\n"
            ),
            t_global,
            iter,
            Est$Q,
            best_Q,
            ifelse(
              is.na(best_iter),
              "NA",
              best_iter
            ),
            sigma,
            relchg,
            dec_count,
            L
          )
        )

      } else {

        cat(
          sprintf(
            paste0(
              "t=%d iter=%d sigma=%.6f ",
              "relchg=%.3e |S_t|=%d\n"
            ),
            t_global,
            iter,
            sigma,
            relchg,
            L
          )
        )
      }
    }

    # (A) parameter-stability stopping rule
    if (relchg < epsilon) {

      stop_reason <- "parameter_convergence"

      if (compute_Q) {
        Q_path <- Q_path[seq_len(iter)]
      }

      break
    }

    # (B) original Q-decrease safeguard
    if (
      compute_Q &&
        dec_count >= q_patience
    ) {

      stop_reason <- "q_decrease"

      Q_path <- Q_path[seq_len(iter)]

      if (verbose) {
        cat(
          sprintf(
            paste0(
              "t=%d early-stop: Q decreased ",
              "%d consecutive times (iter=%d). ",
              "Return best-Q parameters at iter=%d.\n"
            ),
            t_global,
            dec_count,
            iter,
            best_iter
          )
        )
      }

      break
    }
  }

  if (compute_Q) {

    list(
      theta = best_theta,
      Z = best_Z,
      sigma = best_sigma,
      Q_path = Q_path,
      best_Q = best_Q,
      best_iter = best_iter,
      iterations = iter,
      stop_reason = stop_reason,
      neigh_idx = tw$idx,
      neigh_w = tw$w
    )

  } else {

    list(
      theta = theta,
      Z = Z,
      sigma = sigma,
      Q_path = Q_path,
      best_Q = NULL,
      best_iter = NULL,
      iterations = iter,
      stop_reason = stop_reason,
      neigh_idx = tw$idx,
      neigh_w = tw$w
    )
  }
}


# ============================================================
# LA-EM over the full time path
# ============================================================
fit_laem_path <- function(
    Y,
    c,
    h,
    kernel = kernel_epanechnikov,
    self_loop = FALSE,
    smooth = TRUE,
    leave_one_time_out = FALSE,
    max_iter = 1000,
    Q = 10,
    epsilon = 1e-5,
    damp = 0.8,
    verbose = FALSE,
    compute_Q = TRUE,
    omega_newton_maxit = 20,
    omega_newton_tol = 1e-5,
    initial_values,
    n_cores = NULL
) {

  .medcsbm_require("statmod")

  if (length(dim(Y)) != 3L) {
    stop(
      "Y must be an M x E x T edge-only array.",
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
  E <- dim(Y)[2]
  nT <- dim(Y)[3]
  G <- length(unique(c))

  if (
    smooth &&
      (!is.finite(h) || h <= 0)
  ) {
    stop(
      "For smoothed LA-EM, h must be a positive finite number.",
      call. = FALSE
    )
  }

  eg <- prepare_edge_index(
    n,
    self_loop = self_loop
  )

  if (eg$E != E) {
    stop(
      sprintf(
        paste0(
          "Y has E=%d edge columns, but n=%d ",
          "implies E=%d upper-triangular edges."
        ),
        E,
        n,
        eg$E
      ),
      call. = FALSE
    )
  }

  if (!smooth) {
    leave_one_time_out <- FALSE
  }

  if (
    missing(initial_values) ||
      is.null(initial_values)
  ) {
    stop(
      "initial_values must be supplied to fit_laem_path().",
      call. = FALSE
    )
  }

  if (
    !identical(
      dim(initial_values$theta),
      c(n, nT)
    )
  ) {
    stop(
      "initial_values$theta must have dimension n x T.",
      call. = FALSE
    )
  }

  if (
    !identical(
      dim(initial_values$Z),
      c(G, G, nT)
    )
  ) {
    stop(
      "initial_values$Z must have dimension G x G x T.",
      call. = FALSE
    )
  }

  if (length(initial_values$sigma) != nT) {
    stop(
      "initial_values$sigma must have length T.",
      call. = FALSE
    )
  }

  tw_list <- build_time_neighborhoods(
    nT = nT,
    h = h,
    K = kernel,
    smooth = smooth,
    leave_one_time_out = leave_one_time_out
  )

  theta_init <- matrix(
    NA_real_,
    nrow = n,
    ncol = nT
  )

  Z_init <- array(
    NA_real_,
    dim = c(G, G, nT)
  )

  sigma_init <- numeric(nT)

  for (tt in seq_len(nT)) {

    theta_init[, tt] <- center_theta_by_community(
      initial_values$theta[, tt],
      c
    )

    Z_init[, , tt] <- symmetrize_Z(
      initial_values$Z[, , tt]
    )

    sigma_init[tt] <- max(
      as.numeric(initial_values$sigma[tt]),
      1e-8
    )
  }

  jobs <- lapply(
    seq_len(nT),
    function(tt) {
      list(
        t = tt,
        tw = tw_list[[tt]],
        theta0 = theta_init[, tt],
        Z0 = Z_init[, , tt],
        sigma0 = sigma_init[tt]
      )
    }
  )

  fit_job <- function(job) {

    tt <- job$t
    tw <- job$tw

    Y_sub <- Y[
      ,
      ,
      tw$idx,
      drop = FALSE
    ]

    fit_t <- fit_laem_time(
      Y_sub = Y_sub,
      c = c,
      t_global = tt,
      tw = tw,
      eg = eg,
      theta = job$theta0,
      Z = job$Z0,
      sigma = job$sigma0,
      max_iter = max_iter,
      Q = Q,
      epsilon = epsilon,
      damp = damp,
      verbose = verbose,
      compute_Q = compute_Q,
      omega_newton_maxit = omega_newton_maxit,
      omega_newton_tol = omega_newton_tol
    )

    list(
      t = tt,
      fit = fit_t
    )
  }

  n_workers <- resolve_n_cores(
    n_cores,
    n_tasks = nT
  )

  if (n_workers <= 1L) {

    res_list <- lapply(
      jobs,
      fit_job
    )

  } else {

    cl <- parallel::makePSOCKcluster(
      n_workers
    )

    on.exit(
      try(
        parallel::stopCluster(cl),
        silent = TRUE
      ),
      add = TRUE
    )

    local_vars <- c(
      "Y",
      "c",
      "eg",
      "max_iter",
      "Q",
      "epsilon",
      "damp",
      "verbose",
      "compute_Q",
      "omega_newton_maxit",
      "omega_newton_tol"
    )

    parallel::clusterExport(
      cl,
      local_vars,
      envir = environment()
    )

    helper_vars <- c(
      ".medcsbm_require",
      "logistic_prob",
      "log1pexp_stable",
      "center_theta_by_community",
      "symmetrize_Z",
      "compute_edge_linear_predictor",
      "make_block_index",
      "posterior_mode_variance",
      "score_from_residuals",
      "diagonal_information_from_weights",
      "theta_constrained_step",
      "laem_e_step",
      "fit_laem_time"
    )

    parallel::clusterExport(
      cl,
      helper_vars,
      envir = environment(fit_laem_time)
    )

    res_list <- parallel::parLapply(
      cl,
      jobs,
      fit_job
    )

    parallel::stopCluster(cl)
  }

  res_list <- res_list[
    order(
      vapply(
        res_list,
        function(x) x$t,
        integer(1L)
      )
    )
  ]

  tvtheta <- matrix(
    NA_real_,
    nrow = n,
    ncol = nT
  )

  tvZ <- array(
    NA_real_,
    dim = c(G, G, nT)
  )

  tvSigma <- numeric(nT)

  diagnostics <- vector(
    "list",
    nT
  )

  for (ii in seq_along(res_list)) {

    tt <- res_list[[ii]]$t

    tvtheta[, tt] <-
      res_list[[ii]]$fit$theta

    tvZ[, , tt] <-
      res_list[[ii]]$fit$Z

    tvSigma[tt] <-
      res_list[[ii]]$fit$sigma

    diagnostics[[tt]] <-
      res_list[[ii]]$fit[
        c(
          "Q_path",
          "best_Q",
          "best_iter",
          "iterations",
          "stop_reason",
          "neigh_idx",
          "neigh_w"
        )
      ]
  }

  list(
    theta = tvtheta,
    Z = tvZ,
    sigma = tvSigma,
    sigma2 = tvSigma^2,
    diagnostics = diagnostics,
    settings = list(
      smooth = smooth,
      leave_one_time_out = leave_one_time_out,
      Q = Q,
      damp = damp,
      epsilon = epsilon,
      max_iter = max_iter,
      h = h,
      self_loop = self_loop,
      n_cores = n_workers,
      compute_Q = compute_Q,
      q_patience = 5L,
      theta_constraint = "KKT"
    )
  )
}
