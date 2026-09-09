# ============================================================
# inference_core.R
# Pointwise inference for the full-sample LA-EM fit.
#
# Contains:
#   - tangent-space basis under C^T theta = 0 and Z = Z^T
#   - projected score / observed information
#   - sandwich covariance
#   - O(h^2) bias estimate
#   - pointwise Wald confidence intervals
#
# Depends on: base_functions.R, laem_core.R
# ============================================================

# ---- Tangent-space sandwich inference ----
safe_inverse <- function(A, ridge = 0){
  A <- as.matrix(A)
  if(ridge > 0) A <- A + diag(ridge, nrow(A))
  out <- tryCatch(solve(A), error = function(e) NULL)
  if(!is.null(out)) return(out)
  
  ee <- eigen(symmetrize_matrix(A), symmetric = TRUE)
  vals <- ee$values
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(vals)))
  keep <- abs(vals) > tol
  if(!any(keep)) return(matrix(0, nrow(A), ncol(A)))
  
  ee$vectors[, keep, drop = FALSE] %*%
    (diag(1 / vals[keep], sum(keep)) %*% t(ee$vectors[, keep, drop = FALSE]))
}

kernel_moment <- function(kernel, x_power, y_power,
                          lower = -1, upper = 1,
                          subdivisions = 1000){
  stats::integrate(
    function(u) u^x_power * kernel(u)^y_power,
    lower = lower, upper = upper,
    subdivisions = subdivisions, rel.tol = 1e-8
  )$value
}

tangent_space_basis <- function(c, tol = 1e-10){
  c <- as.integer(factor(c, levels = sort(unique(c))))
  n <- length(c)
  G <- length(unique(c))
  p <- n + G * G + 1L
  
  idx <- list(
    n = n, G = G, p_xi = p,
    theta = seq_len(n),
    Z = matrix(n + seq_len(G * G), G, G),
    sigma2 = p,
    c = c
  )
  
  Cmat <- model.matrix(~ factor(c, levels = seq_len(G)) - 1)
  A <- matrix(0, G + G * (G - 1L) / 2L, p)
  A[seq_len(G), idx$theta] <- t(Cmat)
  
  row <- G
  if(G >= 2L){
    for(g in seq_len(G - 1L)){
      for(l in (g + 1L):G){
        row <- row + 1L
        A[row, idx$Z[g, l]] <- 1
        A[row, idx$Z[l, g]] <- -1
      }
    }
  }
  
  qrA <- qr(t(A), tol = tol, LAPACK = FALSE)
  rankA <- qrA$rank
  Qfull <- qr.Q(qrA, complete = TRUE)
  P <- if(rankA >= p) matrix(0, p, 0L) else Qfull[, (rankA + 1L):p, drop = FALSE]
  if(rankA == 0L) P <- diag(p)
  
  list(P = P, A = A, rank = rankA, idx = idx)
}

parameter_curve <- function(fit){
  nT <- ncol(fit$theta)
  p <- nrow(fit$theta) + prod(dim(fit$Z)[1:2]) + 1L
  out <- matrix(NA_real_, p, nT)
  for(t in seq_len(nT)){
    out[, t] <- c(as.numeric(fit$theta[, t]), as.numeric(fit$Z[, , t]), fit$sigma[t]^2)
  }
  out
}

build_inference_context <- function(Y, c, t, fit,
                    h = NULL,
                    kernel = kernel_uniform,
                    self_loop = FALSE,
                    smooth = NULL,
                    leave_one_time_out = NULL,
                    kernel_support = c(-1, 1)){
  stopifnot(length(dim(Y)) == 3)
  M <- dim(Y)[1]
  E <- dim(Y)[2]
  nT <- dim(Y)[3]
  
  c <- as.integer(factor(c, levels = sort(unique(c))))
  eg <- prepare_edge_index(length(c), self_loop)
  stopifnot(eg$E == E)
  
  if(is.null(h)) h <- fit$settings$h
  if(is.null(smooth)) smooth <- isTRUE(fit$settings$smooth)
  if(is.null(leave_one_time_out)){
    leave_one_time_out <- isTRUE(fit$settings$leave_one_time_out)
  }
  
  tw <- if(!smooth){
    list(idx = t, w = 1)
  } else if(leave_one_time_out){
    compute_time_weights_loo(t, nT, h, kernel)
  } else {
    compute_time_weights(t, nT, h, kernel)
  }
  
  basis <- tangent_space_basis(c)
  idx <- basis$idx
  xi_hat <- c(as.numeric(fit$theta[, t]), as.numeric(fit$Z[, , t]), fit$sigma[t]^2)
  theta <- xi_hat[idx$theta]
  Z <- matrix(xi_hat[as.vector(idx$Z)], idx$G, idx$G)
  sigma2 <- max(xi_hat[idx$sigma2], .Machine$double.eps)
  
  list(
    t = t, Y = Y, M = M, E = E, nT = nT, c = c, n = length(c),
    eg = eg, tw = tw, basis = basis, P = basis$P, idx = idx,
    q = ncol(basis$P), xi_hat = xi_hat,
    theta_hat = theta, Z_hat = Z, sigma2_hat = sigma2,
    base_eta = compute_edge_linear_predictor(theta, Z, c, eg$i_idx, eg$j_idx),
    h = h, smooth = smooth,
    nu21 = kernel_moment(kernel, 2, 1, kernel_support[1], kernel_support[2]),
    nu02 = kernel_moment(kernel, 0, 2, kernel_support[1], kernel_support[2]),
    kernel_support = kernel_support
  )
}

sum_by_group <- function(x, group, nbins){
  out <- numeric(nbins)
  if(length(x) == 0L) return(out)
  
  rs <- rowsum(as.numeric(x), as.integer(group), reorder = FALSE)
  id <- as.integer(rownames(rs))
  out[id] <- as.numeric(rs)
  out
}

edge_score_vector <- function(r_edge, eg, c, idx){
  score <- numeric(idx$p_xi)
  score[idx$theta] <- sum_by_group(r_edge, eg$i_idx, idx$n) +
    sum_by_group(r_edge, eg$j_idx, idx$n)
  
  z_id <- c[eg$i_idx] + (c[eg$j_idx] - 1L) * idx$G
  score[as.vector(idx$Z)] <- sum_by_group(r_edge, z_id, idx$G * idx$G)
  score
}

fixed_effect_information <- function(w_edge, eg, c, idx){
  H <- matrix(0, idx$p_xi, idx$p_xi)
  z_pos <- idx$Z[cbind(c[eg$i_idx], c[eg$j_idx])]
  
  for(e in seq_along(w_edge)){
    if(eg$i_idx[e] == eg$j_idx[e]){
      pos <- c(idx$theta[eg$i_idx[e]], z_pos[e])
      x <- c(2, 1)
    } else {
      pos <- c(idx$theta[eg$i_idx[e]], idx$theta[eg$j_idx[e]], z_pos[e])
      x <- c(1, 1, 1)
    }
    H[pos, pos] <- H[pos, pos] + w_edge[e] * tcrossprod(x)
  }
  H
}

subject_score_information <- function(y_k, 
                               base_eta, 
                               sigma2, 
                               eg, 
                               c, 
                               idx,
                               Q = 15,
                               omega_newton_maxit = 50,
                               omega_newton_tol = 1e-8){
  sigma2 <- max(as.numeric(sigma2), .Machine$double.eps)
  sigma <- sqrt(sigma2)
  
  pm <- posterior_mode_variance(
    y_k = y_k,
    base_eta = base_eta,
    sigma = sigma,
    maxit = omega_newton_maxit,
    tol = omega_newton_tol
  )
  
  mode <- pm$mean
  v_post <- max(pm$var, 1e-12)
  
  gh <- statmod::gauss.quad(Q, "hermite")
  omega <- mode + sqrt(2 * v_post) * gh$nodes
  
  log_joint <- numeric(Q)
  for(q in seq_len(Q)){
    eta_q <- base_eta + omega[q]
    log_joint[q] <- sum(y_k * eta_q - log1pexp_stable(eta_q)) -
      0.5 * (log(2 * pi * sigma2) + omega[q]^2 / sigma2)
  }
  
  log_post_w <- log(gh$weights) + log_joint + gh$nodes^2
  lambda <- exp(log_post_w - log_sum_exp(log_post_w))
  
  score_nodes <- matrix(0, nrow = Q, ncol = idx$p_xi)
  w_bar <- numeric(length(base_eta))
  neg_hess_sigma_nodes <- numeric(Q)
  
  for(q in seq_len(Q)){
    eta_q <- base_eta + omega[q]
    mu_q <- logistic_prob(eta_q)
    r_q <- y_k - mu_q
    w_q <- mu_q * (1 - mu_q)
    
    score_nodes[q, ] <- edge_score_vector(r_q, eg, c, idx)
    score_nodes[q, idx$sigma2] <-
      -0.5 / sigma2 + 0.5 * omega[q]^2 / sigma2^2
    
    w_bar <- w_bar + lambda[q] * w_q
    
    neg_hess_sigma_nodes[q] <-
      omega[q]^2 / sigma2^3 - 0.5 / sigma2^2
  }
  
  score <- as.vector(crossprod(lambda, score_nodes))
  
  centered <- sweep(score_nodes, 2L, score, "-")
  cov_score <- crossprod(centered * sqrt(lambda),
                         centered * sqrt(lambda))
  
  neg_complete_info <- fixed_effect_information(w_bar, eg, c, idx)
  neg_complete_info[idx$sigma2, idx$sigma2] <-
    sum(lambda * neg_hess_sigma_nodes)
  
  neg_info <- neg_complete_info - cov_score
  neg_info <- 0.5 * (neg_info + t(neg_info))
  
  list(score = score, neg_info = neg_info)
}

time_score_information <- function(ctx,
                                   smooth = ctx$smooth,
                                   Q = 15,
                                   omega_newton_maxit = 50,
                                   omega_newton_tol = 1e-8,
                                   center_scores = TRUE,
                                   compute_meat = TRUE){
  H_R <- matrix(0, ctx$q, ctx$q)
  meat_R <- if(compute_meat && !center_scores) matrix(0, ctx$q, ctx$q) else NULL
  score_mat <- meat_w <- NULL
  if(compute_meat && center_scores){
    nr <- ctx$M * length(ctx$tw$idx)
    score_mat <- matrix(0, nr, ctx$q)
    meat_w <- numeric(nr)
    rr <- 0L
  }
  for(ss in seq_along(ctx$tw$idx)){
    s <- ctx$tw$idx[ss]
    ws <- ctx$tw$w[ss]
    Y_s <- ctx$Y[, , s, drop = FALSE][, , 1]
    for(m in seq_len(ctx$M)){
      obs <- subject_score_information(
        y_k = Y_s[m, ], base_eta = ctx$base_eta,
        sigma2 = ctx$sigma2_hat, eg = ctx$eg, c = ctx$c, idx = ctx$idx,
        Q = Q, omega_newton_maxit = omega_newton_maxit,
        omega_newton_tol = omega_newton_tol
      )
      score_R <- as.vector(crossprod(ctx$P, obs$score))
      H_R <- H_R + ws * crossprod(ctx$P, obs$neg_info %*% ctx$P)
      if(compute_meat && center_scores){
        rr <- rr + 1L
        score_mat[rr, ] <- score_R
        meat_w[rr] <- ws^2
      } else if(compute_meat){
        meat_R <- meat_R + ws^2 * tcrossprod(score_R)
      }
    }
  }
  J_hat <- H_R / ctx$M
  out <- list(H_R = H_R, J_hat = J_hat,
              score_mat = score_mat, meat_weights = meat_w)
  if(compute_meat){
    if(center_scores){
      score_bar <- colSums(score_mat * meat_w) / sum(meat_w)
      score_centered <- sweep(score_mat, 2L, score_bar, "-")
      meat_R <- crossprod(score_centered * sqrt(meat_w))
      out$score_centered <- score_centered
      out$score_bar <- score_bar
    }
    out$meat_R <- meat_R
    out$Upsilon_hat <- if(smooth) (ctx$nT * ctx$h / ctx$nu02) * (meat_R / ctx$M) else meat_R / ctx$M
    out$Upsilon_type <- if(smooth) "kernel_weighted_pointwise_rescaled" else "pointwise"
  }
  out
}

estimate_variance_time <- function(Y, c, t, fit,
                          h = NULL,
                          kernel = kernel_uniform,
                          self_loop = FALSE,
                          smooth = TRUE,
                          Q = 15,
                          kernel_support = c(-1, 1),
                          omega_newton_maxit = 50,
                          omega_newton_tol = 1e-8,
                          ridge = 1e-8,
                          center_scores = TRUE,
                          return_vcov = TRUE){
  ctx <- build_inference_context(Y, c, t, fit, h, kernel, self_loop, smooth, FALSE, kernel_support)
  info <- time_score_information(
    ctx, smooth, Q, omega_newton_maxit, omega_newton_tol,
    center_scores = center_scores, compute_meat = TRUE
  )
  
  J_inv <- safe_inverse(info$J_hat, ridge)
  core <- ctx$P %*% J_inv %*% info$Upsilon_hat %*% J_inv %*% t(ctx$P)
  Sigma_xi <- symmetrize_matrix(if(smooth) ctx$nu02 * core else core)
  scale <- if(smooth) ctx$M * ctx$nT * ctx$h else ctx$M
  vcov_xi <- symmetrize_matrix(Sigma_xi / scale)
  var_xi <- pmax(diag(vcov_xi), 0)
  se_xi <- sqrt(var_xi)
  
  out <- list(
    t = t, smooth = smooth, neigh_idx = ctx$tw$idx, neigh_w = ctx$tw$w,
    xi_hat = ctx$xi_hat, H_R = info$H_R, meat_R = info$meat_R,
    J_hat = info$J_hat, Upsilon_hat = info$Upsilon_hat,
    Sigma_xi = Sigma_xi, var_xi = var_xi, se_xi = se_xi,
    theta_var = var_xi[ctx$idx$theta],
    Z_var = matrix(var_xi[as.vector(ctx$idx$Z)], ctx$idx$G, ctx$idx$G),
    sigma2_var = var_xi[ctx$idx$sigma2],
    theta_se = se_xi[ctx$idx$theta],
    Z_se = matrix(se_xi[as.vector(ctx$idx$Z)], ctx$idx$G, ctx$idx$G),
    sigma2_se = se_xi[ctx$idx$sigma2],
    P = ctx$P, constraints = ctx$basis$A, nu02 = ctx$nu02,
    idx = ctx$idx, score_mat = info$score_mat,
    meat_weights = info$meat_weights,
    score_centered = info$score_centered,
    settings = list(
      h = ctx$h, Q = Q, ridge = ridge,
      kernel_support = kernel_support, smooth = smooth,
      center_scores = center_scores,
      J = if(smooth) "kernel_weighted" else "pointwise",
      Upsilon = info$Upsilon_type
    )
  )
  
  if(return_vcov) out$vcov_xi <- vcov_xi
  out
}

estimate_second_derivative <- function(fit, t, kernel = kernel_uniform, dh = NULL){
  if(is.null(dh)) dh <- fit$settings$h
  X <- parameter_curve(fit)
  nT <- ncol(X)
  if(nT < 3L) stop("Need at least 3 time points.")
  
  u <- seq_len(nT) / nT
  x <- u - u[t]
  w <- kernel(x / dh)
  w[!is.finite(w) | w < 0] <- 0
  keep <- which(w > 0)
  if(length(keep) < 3L){
    keep <- order(abs(x))[seq_len(3L)]
    w <- rep(0, nT)
    w[keep] <- 1
  }
  
  D <- cbind(1, x[keep], 0.5 * x[keep]^2)
  Y <- t(X[, keep, drop = FALSE])
  sw <- sqrt(w[keep])
  as.vector(qr.solve(D * sw, Y * sw)[3L, ])
}

estimate_smoothing_bias <- function(P, fit, t = NULL, h = NULL,
                     kernel = kernel_uniform,
                     support = c(-1, 1)){
  if(is.list(P) && !is.null(P$P)){
    v <- P
    P <- v$P
    if(is.null(t)) t <- v$t
    if(is.null(h)) h <- v$settings$h
  }
  if(is.null(t) || is.null(h)) stop("Need t and h.")
  
  x2 <- estimate_second_derivative(fit, t, kernel, h)
  nu21 <- kernel_moment(kernel, 2L, 1L, support[1], support[2])
  b <- as.vector(0.5 * nu21 * h^2 * P %*% crossprod(P, x2))
  list(bias = b, estimate_second_derivative = x2, nu21 = nu21)
}

construct_confidence_interval <- function(fit, t, bias, se, level = 0.95) {
  est <- c(
    as.numeric(fit$theta[, t]),
    as.numeric(fit$Z[, , t]),
    fit$sigma[t]^2
  )
  if (is.list(bias)) bias <- bias$bias
  if (is.list(se)) se <- se$se_xi
  stopifnot(length(est) == length(bias), length(est) == length(se))

  half_width <- stats::qnorm(1 - (1 - level) / 2) * se
  estimate_bc <- est - bias

  cbind(
    estimate = est,
    bias = bias,
    se = se,
    estimate_bc = estimate_bc,
    half_width = half_width,
    ci_width = 2 * half_width,
    lower = estimate_bc - half_width,
    upper = estimate_bc + half_width
  )
}


