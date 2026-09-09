# ============================================================
# auxiliary_model.R
# Community-free auxiliary fit used by RJSE
# ============================================================
#
# Pointwise Bernoulli working model:
#
#   logit P{Y_m,ij(u^t)=1}
#     = theta_i(u^t) + theta_j(u^t) + omega_m(u^t),
#
# with sum_m omega_m(u^t) = 0.
#
# The fitted probabilities are used only to construct Pearson
# residuals for community recovery.  These auxiliary parameters
# are not the mixed-effects DCSBM parameters.
# ============================================================


# ------------------------------------------------------------
# Infer n from the number of edge columns
# ------------------------------------------------------------
rjse_infer_n_from_edges <- function(E, self_loop = FALSE) {

  if (self_loop) {
    n <- (sqrt(1 + 8 * E) - 1) / 2
  } else {
    n <- (1 + sqrt(1 + 8 * E)) / 2
  }

  n_int <- as.integer(round(n))

  if (abs(n - n_int) > 1e-8) {
    stop(
      "Cannot infer the number of nodes from the edge dimension.",
      call. = FALSE
    )
  }

  n_int
}


# ------------------------------------------------------------
# Identifiability for subject effects
# ------------------------------------------------------------
rjse_center_omega <- function(theta, omega) {

  omega_bar <- mean(omega)

  list(
    theta = theta + 0.5 * omega_bar,
    omega = omega - omega_bar
  )
}


# ------------------------------------------------------------
# Pointwise auxiliary initialization
# ------------------------------------------------------------
rjse_aux_initialization <- function(Y_t, n, pseudo = 0.5) {

  stopifnot(is.matrix(Y_t))

  M <- nrow(Y_t)
  E <- ncol(Y_t)

  p_hat <- (sum(Y_t) + pseudo) /
    (M * E + 2 * pseudo)

  list(
    theta = rep(0.5 * stats::qlogis(p_hat), n),
    omega = numeric(M)
  )
}


# ------------------------------------------------------------
# KKT step for sum_m omega_m = 0
# ------------------------------------------------------------
rjse_omega_kkt <- function(
    S_omega,
    H_pos_omega,
    damp = 1.0,
    eps = 1e-15
) {

  H_pos_omega <- pmax(H_pos_omega, eps)
  invH <- 1 / H_pos_omega
  b <- damp * invH * S_omega

  alpha <- sum(b) / pmax(sum(invH), eps)

  b - invH * alpha
}


rjse_cap_step <- function(step, max_step = Inf) {
  if (!is.finite(max_step)) return(step)
  pmax(pmin(step, max_step), -max_step)
}


# ------------------------------------------------------------
# Score / diagonal information for the auxiliary model
# ------------------------------------------------------------
rjse_aux_stats <- function(
    y_list,
    w_time,
    theta,
    omega,
    i_idx,
    j_idx,
    n,
    edge_chunk_size = 50000L,
    compute_loglik = TRUE
) {

  M <- nrow(y_list[[1]])
  E <- length(i_idx)

  if (!is.finite(edge_chunk_size) || edge_chunk_size <= 0) {
    edge_chunk_size <- E
  }

  edge_chunk_size <- as.integer(min(E, edge_chunk_size))

  g_theta <- numeric(n)
  H_theta <- numeric(n)
  g_omega <- numeric(M)
  H_omega <- numeric(M)
  ll <- 0

  for (s in seq_along(y_list)) {

    y_mat_s <- y_list[[s]]
    ws <- w_time[s]

    for (start in seq.int(1L, E, by = edge_chunk_size)) {

      end <- min(E, start + edge_chunk_size - 1L)
      idx <- start:end

      y_chunk <- y_mat_s[, idx, drop = FALSE]

      base_eta <- theta[i_idx[idx]] + theta[j_idx[idx]]

      eta <- matrix(
        base_eta,
        nrow = M,
        ncol = length(idx),
        byrow = TRUE
      )

      eta <- sweep(eta, 1L, omega, "+")

      p <- stats::plogis(eta)
      R <- y_chunk - p
      W <- p * (1 - p)

      r_edge <- ws * colSums(R)
      w_edge <- ws * colSums(W)

      tmp_i <- rowsum(r_edge, i_idx[idx], reorder = FALSE)
      tmp_j <- rowsum(r_edge, j_idx[idx], reorder = FALSE)

      ii <- as.integer(rownames(tmp_i))
      jj <- as.integer(rownames(tmp_j))

      g_theta[ii] <- g_theta[ii] + as.numeric(tmp_i)
      g_theta[jj] <- g_theta[jj] + as.numeric(tmp_j)

      tmp_i <- rowsum(w_edge, i_idx[idx], reorder = FALSE)
      tmp_j <- rowsum(w_edge, j_idx[idx], reorder = FALSE)

      ii <- as.integer(rownames(tmp_i))
      jj <- as.integer(rownames(tmp_j))

      H_theta[ii] <- H_theta[ii] + as.numeric(tmp_i)
      H_theta[jj] <- H_theta[jj] + as.numeric(tmp_j)

      same <- i_idx[idx] == j_idx[idx]

      if (any(same)) {
        tmp_same <- rowsum(
          2 * w_edge[same],
          i_idx[idx][same],
          reorder = FALSE
        )

        ss <- as.integer(rownames(tmp_same))
        H_theta[ss] <- H_theta[ss] + as.numeric(tmp_same)
      }

      g_omega <- g_omega + ws * rowSums(R)
      H_omega <- H_omega + ws * rowSums(W)

      if (compute_loglik) {
        ll <- ll + ws * sum(
          y_chunk * eta - log1pexp_stable(eta)
        )
      }
    }
  }

  list(
    g_theta = g_theta,
    H_theta = H_theta,
    g_omega = g_omega,
    H_omega = H_omega,
    logLik = ll
  )
}


# ------------------------------------------------------------
# Pearson residuals at one time point
# ------------------------------------------------------------
rjse_pearson_residual_time <- function(
    Y_t,
    theta,
    omega,
    eg,
    eps = 1e-6,
    edge_chunk_size = 50000L
) {

  stopifnot(is.matrix(Y_t))

  M <- nrow(Y_t)
  E <- ncol(Y_t)

  i_idx <- eg$i_idx
  j_idx <- eg$j_idx

  stopifnot(
    length(i_idx) == E,
    length(j_idx) == E
  )

  if (!is.finite(edge_chunk_size) || edge_chunk_size <= 0) {
    edge_chunk_size <- E
  }

  edge_chunk_size <- as.integer(min(E, edge_chunk_size))

  P0 <- matrix(NA_real_, nrow = M, ncol = E)
  R <- matrix(NA_real_, nrow = M, ncol = E)

  for (start in seq.int(1L, E, by = edge_chunk_size)) {

    end <- min(E, start + edge_chunk_size - 1L)
    idx <- start:end

    base_eta <- theta[i_idx[idx]] + theta[j_idx[idx]]

    eta <- matrix(
      base_eta,
      nrow = M,
      ncol = length(idx),
      byrow = TRUE
    )

    eta <- sweep(eta, 1L, omega, "+")

    p <- pmin(
      pmax(stats::plogis(eta), eps),
      1 - eps
    )

    P0[, idx] <- p

    R[, idx] <-
      (Y_t[, idx, drop = FALSE] - p) /
      sqrt(p * (1 - p))
  }

  list(
    R = R,
    P0 = P0
  )
}


# ============================================================
# Auxiliary fit at one time point
# ============================================================
rjse_fit_auxiliary_time <- function(
    Y_sub,
    t_global,
    tw,
    eg,
    theta = NULL,
    omega = NULL,
    max_iter = 300,
    epsilon = 1e-5,
    damp = 0.8,
    verbose = FALSE,
    edge_chunk_size = 50000L,
    max_step = 5,
    ll_patience = 5L,
    ll_decrease_tol = 1e-10
) {

  stopifnot(length(dim(Y_sub)) == 3L)

  M <- dim(Y_sub)[1]
  E <- dim(Y_sub)[2]
  L <- dim(Y_sub)[3]

  stopifnot(
    length(tw$idx) == L,
    length(tw$w) == L
  )

  i_idx <- eg$i_idx
  j_idx <- eg$j_idx

  stopifnot(
    length(i_idx) == E,
    length(j_idx) == E
  )

  n <- max(i_idx, j_idx)

  y_list <- lapply(
    seq_len(L),
    function(s) {
      Y_sub[, , s, drop = FALSE][, , 1]
    }
  )

  w_s <- tw$w

  if (is.null(theta) || is.null(omega)) {

    y_bar <- Reduce(
      `+`,
      Map(
        function(y, w) w * y,
        y_list,
        w_s
      )
    )

    init <- rjse_aux_initialization(
      y_bar,
      n = n,
      pseudo = 0.5
    )

    if (is.null(theta)) theta <- init$theta
    if (is.null(omega)) omega <- init$omega
  }

  centered <- rjse_center_omega(theta, omega)
  theta <- centered$theta
  omega <- centered$omega

  best_ll <- -Inf
  best_iter <- NA_integer_
  best_theta <- theta
  best_omega <- omega

  dec_count <- 0L
  ll_path <- rep(NA_real_, max_iter)

  for (iter in seq_len(max_iter)) {

    st <- rjse_aux_stats(
      y_list = y_list,
      w_time = w_s,
      theta = theta,
      omega = omega,
      i_idx = i_idx,
      j_idx = j_idx,
      n = n,
      edge_chunk_size = edge_chunk_size,
      compute_loglik = FALSE
    )

    theta_step <-
      damp * st$g_theta /
      pmax(st$H_theta, 1e-15)

    omega_step <- rjse_omega_kkt(
      S_omega = st$g_omega,
      H_pos_omega = st$H_omega,
      damp = damp
    )

    theta_new <- theta + rjse_cap_step(
      theta_step,
      max_step = max_step
    )

    omega_new <- omega + rjse_cap_step(
      omega_step,
      max_step = max_step
    )

    centered <- rjse_center_omega(
      theta_new,
      omega_new
    )

    theta_new <- centered$theta
    omega_new <- centered$omega

    relchg <- max(
      abs(
        c(
          theta_new - theta,
          omega_new - omega
        )
      )
    ) /
      (
        1 +
          max(
            abs(
              c(theta, omega)
            )
          )
      )

    theta <- theta_new
    omega <- omega_new

    st_ll <- rjse_aux_stats(
      y_list = y_list,
      w_time = w_s,
      theta = theta,
      omega = omega,
      i_idx = i_idx,
      j_idx = j_idx,
      n = n,
      edge_chunk_size = edge_chunk_size,
      compute_loglik = TRUE
    )

    ll_now <- st_ll$logLik
    ll_path[iter] <- ll_now

    if (is.finite(ll_now) && ll_now > best_ll) {
      best_ll <- ll_now
      best_iter <- iter
      best_theta <- theta
      best_omega <- omega
    }

    if (iter == 1L) {
      dec_count <- 0L
    } else {
      ll_prev <- ll_path[iter - 1L]

      if (
        is.finite(ll_prev) &&
        ll_now < ll_prev - ll_decrease_tol
      ) {
        dec_count <- dec_count + 1L
      } else {
        dec_count <- 0L
      }
    }

    if (verbose) {
      cat(
        sprintf(
          paste0(
            "t=%d iter=%d ll=%.6f best=%.6f(best@%s) ",
            "relchg=%.3e dec=%d |S_t|=%d\n"
          ),
          t_global,
          iter,
          ll_now,
          best_ll,
          ifelse(is.na(best_iter), "NA", best_iter),
          relchg,
          dec_count,
          L
        )
      )
    }

    if (relchg < epsilon) {
      ll_path <- ll_path[seq_len(iter)]
      break
    }

    if (dec_count >= ll_patience) {
      if (verbose) {
        cat(
          sprintf(
            paste0(
              "t=%d early-stop: log-likelihood decreased %d ",
              "consecutive times (iter=%d). Return best-ll ",
              "parameters at iter=%d.\n"
            ),
            t_global,
            dec_count,
            iter,
            best_iter
          )
        )
      }

      ll_path <- ll_path[seq_len(iter)]
      break
    }
  }

  list(
    theta = best_theta,
    omega = best_omega,
    ll_path = ll_path,
    best_ll = best_ll,
    best_iter = best_iter,
    neigh_idx = tw$idx,
    neigh_w = tw$w
  )
}


# ============================================================
# Full pointwise auxiliary fit used by estimate_membership()
# ============================================================
fit_rjse_auxiliary_model <- function(
    Y,
    n = NULL,
    self_loop = FALSE,
    max_iter = 300,
    epsilon = 1e-5,
    damp = 0.8,
    warmstart_max_iter = 30,
    warmstart_epsilon = 1e-3,
    edge_chunk_size = 50000L,
    max_step = 5,
    ll_patience = 5L,
    ll_decrease_tol = 1e-10,
    residual_eps = 1e-6,
    n_cores = NULL,
    verbose = FALSE
) {

  if (length(dim(Y)) != 3L) {
    stop(
      "Y must be an M x E x T edge-only array.",
      call. = FALSE
    )
  }

  if (anyNA(Y) || !all(Y == 0 | Y == 1)) {
    stop(
      "The Bernoulli RJSE auxiliary model requires binary Y.",
      call. = FALSE
    )
  }

  M <- dim(Y)[1]
  E <- dim(Y)[2]
  nT <- dim(Y)[3]

  if (is.null(n)) {
    n <- rjse_infer_n_from_edges(
      E,
      self_loop = self_loop
    )
  }

  eg <- prepare_edge_index(
    n,
    self_loop = self_loop
  )

  if (eg$E != E) {
    stop(
      "Y and n imply different numbers of edges.",
      call. = FALSE
    )
  }

  # The supplied RJSE detection code fits the auxiliary model
  # pointwise, so each target time uses only its own layer.
  tw_list <- lapply(
    seq_len(nT),
    function(tt) {
      list(
        idx = tt,
        w = 1
      )
    }
  )

  # ----------------------------------------------------------
  # Stage 0: short pointwise fit used as the starting value
  # ----------------------------------------------------------
  theta_init <- matrix(NA_real_, n, nT)
  omega_init <- matrix(NA_real_, M, nT)

  for (tt in seq_len(nT)) {

    Y_t <- Y[, , tt, drop = FALSE][, , 1]

    init <- rjse_aux_initialization(
      Y_t,
      n = n,
      pseudo = 0.5
    )

    fit0 <- rjse_fit_auxiliary_time(
      Y_sub = Y[, , tt, drop = FALSE],
      t_global = tt,
      tw = tw_list[[tt]],
      eg = eg,
      theta = init$theta,
      omega = init$omega,
      max_iter = warmstart_max_iter,
      epsilon = warmstart_epsilon,
      damp = damp,
      verbose = FALSE,
      edge_chunk_size = edge_chunk_size,
      max_step = max_step,
      ll_patience = ll_patience,
      ll_decrease_tol = ll_decrease_tol
    )

    theta_init[, tt] <- fit0$theta
    omega_init[, tt] <- fit0$omega
  }

  jobs <- lapply(
    seq_len(nT),
    function(tt) {
      list(
        t = tt,
        theta0 = theta_init[, tt],
        omega0 = omega_init[, tt]
      )
    }
  )

  run_job <- function(job) {

    tt <- job$t
    tw <- tw_list[[tt]]

    fit <- rjse_fit_auxiliary_time(
      Y_sub = Y[, , tt, drop = FALSE],
      t_global = tt,
      tw = tw,
      eg = eg,
      theta = job$theta0,
      omega = job$omega0,
      max_iter = max_iter,
      epsilon = epsilon,
      damp = damp,
      verbose = verbose,
      edge_chunk_size = edge_chunk_size,
      max_step = max_step,
      ll_patience = ll_patience,
      ll_decrease_tol = ll_decrease_tol
    )

    rp <- rjse_pearson_residual_time(
      Y_t = Y[, , tt, drop = FALSE][, , 1],
      theta = fit$theta,
      omega = fit$omega,
      eg = eg,
      eps = residual_eps,
      edge_chunk_size = edge_chunk_size
    )

    list(
      t = tt,
      fit = fit,
      residual = rp
    )
  }

  n_workers <- resolve_n_cores(
    n_cores,
    n_tasks = nT
  )

  if (n_workers <= 1L) {

    res_list <- lapply(
      jobs,
      run_job
    )

  } else {

    cl <- parallel::makePSOCKcluster(n_workers)

    on.exit(
      try(
        parallel::stopCluster(cl),
        silent = TRUE
      ),
      add = TRUE
    )

    local_vars <- c(
      "Y", "tw_list", "eg", "max_iter", "epsilon", "damp",
      "verbose", "edge_chunk_size", "max_step", "ll_patience",
      "ll_decrease_tol", "residual_eps"
    )

    parallel::clusterExport(
      cl,
      local_vars,
      envir = environment()
    )

    helper_vars <- c(
      "log1pexp_stable",
      "rjse_center_omega",
      "rjse_aux_initialization",
      "rjse_omega_kkt",
      "rjse_cap_step",
      "rjse_aux_stats",
      "rjse_pearson_residual_time",
      "rjse_fit_auxiliary_time"
    )

    parallel::clusterExport(
      cl,
      helper_vars,
      envir = environment(rjse_fit_auxiliary_time)
    )

    res_list <- parallel::parLapply(
      cl,
      jobs,
      run_job
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

  tvtheta <- matrix(NA_real_, n, nT)
  tvomega <- matrix(NA_real_, M, nT)

  R_aux <- array(
    NA_real_,
    dim = c(M, E, nT)
  )

  P0_aux <- array(
    NA_real_,
    dim = c(M, E, nT)
  )

  ll_paths <- vector("list", nT)
  best_ll <- rep(NA_real_, nT)
  best_it <- rep(NA_integer_, nT)

  for (ii in seq_along(res_list)) {

    tt <- res_list[[ii]]$t

    tvtheta[, tt] <- res_list[[ii]]$fit$theta
    tvomega[, tt] <- res_list[[ii]]$fit$omega

    R_aux[, , tt] <- res_list[[ii]]$residual$R
    P0_aux[, , tt] <- res_list[[ii]]$residual$P0

    ll_paths[[tt]] <- res_list[[ii]]$fit$ll_path
    best_ll[tt] <- res_list[[ii]]$fit$best_ll
    best_it[tt] <- res_list[[ii]]$fit$best_iter
  }

  list(
    theta = tvtheta,
    omega = tvomega,
    R = R_aux,
    P0 = P0_aux,
    logLik_paths = ll_paths,
    best_logLik = best_ll,
    best_iter = best_it,
    edge_index = eg,
    settings = list(
      model = "community_free_bernoulli_auxiliary",
      M = M,
      n = n,
      E = E,
      T = nT,
      max_iter = max_iter,
      epsilon = epsilon,
      damp = damp,
      residual_eps = residual_eps,
      n_cores = n_workers,
      pointwise = TRUE
    )
  )
}
