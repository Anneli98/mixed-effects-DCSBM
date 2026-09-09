# ============================================================
# marginal_likelihood.R
# Leave-one-time-out bandwidth selection.
#
# For each candidate h and target t:
#   1) exclude the complete time layer t from initialization and LA-EM;
#   2) estimate xi_h^{(-t)}(u^t);
#   3) evaluate the held-out marginal log-likelihood at time t;
#   4) average over subjects and time points.
#
# Depends on: base_functions.R, laem_core.R
# ============================================================

# ---- Predictive marginal likelihood for LOTO CV ----
marginal_loglik_subject <- function(y_k, base_eta, sigma, nAGQ = 10) {
  pm <- posterior_mode_variance(y_k, base_eta, sigma)
  m  <- pm$mean
  v  <- pm$var
  
  gh  <- statmod::gauss.quad(nAGQ, "hermite")
  s   <- gh$nodes
  w   <- gh$weights
  
  sq2v <- sqrt(2 * v)
  omega_q <- m + sq2v * s
  
  val_q <- sapply(seq_along(omega_q), function(q){
    eta_q <- base_eta + omega_q[q]
    ll_edges <- sum(y_k * eta_q - log1pexp_stable(eta_q))
    ll_prior <- -0.5 * ( (omega_q[q]^2) / (sigma^2) + log(2*pi*sigma^2) )
    ll_edges + ll_prior
  })
  
  # note: hermite rule -> integral factor handling
  logI <- log(sq2v) + log_sum_exp( log(w) + val_q + s^2 )
  as.numeric(logI)
}

heldout_marginal_loglik <- function(Y, c, fit, nAGQ = 10, sigma_floor = 1e-8) {
  .medcsbm_require("statmod")

  c <- as.integer(factor(c, levels = sort(unique(c))))
  n <- length(c)
  nT <- dim(Y)[3]
  M <- dim(Y)[1]
  eg <- prepare_edge_index(n, self_loop = isTRUE(fit$settings$self_loop))
  stopifnot(eg$E == dim(Y)[2])

  ll_time <- numeric(nT)
  ll_subject_time <- matrix(NA_real_, M, nT)

  for (t in seq_len(nT)) {
    base_eta <- compute_edge_linear_predictor(
      fit$theta[, t], fit$Z[, , t], c, eg$i_idx, eg$j_idx
    )
    sigma_t <- max(fit$sigma[t], sigma_floor)

    for (m in seq_len(M)) {
      ll_subject_time[m, t] <- marginal_loglik_subject(
        y_k = Y[m, , t],
        base_eta = base_eta,
        sigma = sigma_t,
        nAGQ = nAGQ
      )
    }
    ll_time[t] <- sum(ll_subject_time[, t])
  }

  list(
    total = sum(ll_time),
    mean = sum(ll_time) / (M * nT),
    by_time = ll_time / M,
    by_subject_time = ll_subject_time
  )
}


