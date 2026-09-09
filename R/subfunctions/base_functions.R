# ============================================================
# base_functions.R
# Shared utilities for the time-varying mixed-effects DCSBM.
#
# Contains only functionality used by multiple algorithm modules:
#   - kernels / time weights / edge indexing
#   - identifiability helpers
#   - posterior mode + curvature for the random intercept
#   - tidy parameter metadata
#
# Debug note: source this file first when testing any module alone.
# ============================================================

.medcsbm_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it first.", pkg), call. = FALSE)
  }
}

logistic_prob <- function(x) stats::plogis(x)

log1pexp_stable <- function(z) {
  out <- numeric(length(z))
  pos <- z > 0
  out[pos] <- z[pos] + log1p(exp(-z[pos]))
  out[!pos] <- log1p(exp(z[!pos]))
  dim(out) <- dim(z)
  out
}

kernel_uniform   <- function(u) 0.5 * (abs(u) <= 1)
kernel_epanechnikov <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1)
kernel_triangular    <- function(u) (1 - abs(u)) * (abs(u) <= 1)
kernel_gaussian  <- function(u) stats::dnorm(u)

get_kernel <- function(kernel = "epanechnikov") {
  if (is.function(kernel)) {
    return(list(fun = kernel, name = "custom", support = c(-1, 1)))
  }
  key <- tolower(as.character(kernel)[1])
  if (key %in% c("epanechnikov", "epanech", "epa")) {
    return(list(fun = kernel_epanechnikov, name = "epanechnikov", support = c(-1, 1)))
  }
  if (key %in% c("uniform", "unif")) {
    return(list(fun = kernel_uniform, name = "uniform", support = c(-1, 1)))
  }
  if (key %in% c("triangular", "triangle", "tri")) {
    return(list(fun = kernel_triangular, name = "triangular", support = c(-1, 1)))
  }
  if (key %in% c("gaussian", "normal", "gauss")) {
    return(list(fun = kernel_gaussian, name = "gaussian", support = c(-Inf, Inf)))
  }
  stop("Unknown kernel. Use 'epanechnikov', 'uniform', 'triangular', 'gaussian', or pass a function.",
       call. = FALSE)
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

symmetrize_matrix <- function(A) 0.5 * (A + t(A))


# ---- Data / kernel / constraint helpers ----
prepare_edge_index <- function(n, self_loop = FALSE){
  idx <- upper.tri(matrix(0, n, n), diag = self_loop)
  ij  <- which(idx, arr.ind = TRUE)
  list(
    idx   = idx,
    i_idx = as.integer(ij[,1]),
    j_idx = as.integer(ij[,2]),
    E     = nrow(ij)
  )
}

compute_time_weights <- function(t, nT, h, K) {
  u <- (1:nT)/nT
  w_raw <- K((u - u[t]) / h)
  eps = 1e-15
  w_raw[w_raw < eps] <- 0    
  idx <- which(w_raw > 0)
  if(length(idx) == 0){
    idx <- t; w_raw[t] <- 1; idx <- t
  }
  w <- w_raw[idx]
  w <- w / sum(w)
  list(idx = idx, w = w)
}

compute_time_weights_loo <- function(t, nT, h, K) {
  u <- (1:nT) / nT
  w_raw <- K((u - u[t]) / h)   
  w_raw[t] <- 0               # leave-one-out
  eps = 1e-15
  w_raw[w_raw < eps] <- 0     
  s <- sum(w_raw)
  if (!is.finite(s) || s <= 0) {
    stop(sprintf("LOO weights undefined at t=%d: sum of weights is zero.", t))
  }
  
  idx <- which(w_raw > 0)
  w <- w_raw[idx] / s
  list(idx = idx, w = w)
}

build_time_neighborhoods <- function(nT, h, K,
                         smooth = TRUE,
                         leave_one_time_out = FALSE){
  if(!smooth){
    return(lapply(seq_len(nT), function(t) list(idx = t, w = 1)))
  }
  if(leave_one_time_out){
    return(lapply(seq_len(nT), function(t) compute_time_weights_loo(t, nT, h, K = K)))
  }
  lapply(seq_len(nT), function(t) compute_time_weights(t, nT, h, K = K))
}

center_theta_by_community <- function(theta, c){
  theta - ave(theta, c)
}

symmetrize_Z <- function(Z){
  0.5 * (Z + t(Z))
}

compute_edge_linear_predictor <- function(theta, Z, c, i_idx, j_idx){
  theta[i_idx] + theta[j_idx] + Z[cbind(c[i_idx], c[j_idx])]
}

make_block_index <- function(c, i_idx, j_idx, G){
  c[i_idx] + (c[j_idx] - 1L) * G
}


# ---- Shared posterior approximation ----
posterior_mode_variance <- function(y_k, base_eta, sigma, m0 = 0,
                            maxit = 50, tol = 1e-8,
                            H_floor = 1e-10,
                            max_step = 10,
                            max_backtrack = 20){
  
  logpost <- function(m){
    eta <- base_eta + m
    sum(y_k * eta - log1pexp_stable(eta)) - 0.5 * (m^2) / (sigma^2)
  }
  
  m <- m0
  invs2 <- 1/(sigma^2)
  
  for(it in 1:maxit){
    eta <- base_eta + m
    mu  <- 1/(1+exp(-eta))
    g   <- sum(y_k - mu) - m * invs2
    H   <- -sum(mu*(1-mu)) - invs2  # negative
    if(abs(H) < H_floor) H <- -H_floor
    
    step <- g / H
    step <- max(min(step,  max_step), -max_step)
    m_new <- m - step
    
    lp0 <- logpost(m)
    lp1 <- logpost(m_new)
    bt <- 0
    while(lp1 < lp0 && bt < max_backtrack){
      step <- step / 2
      m_new <- m - step
      lp1 <- logpost(m_new)
      bt <- bt + 1
    }
    
    if(abs(step) < tol){
      m <- m_new
      break
    }
    m <- m_new
  }
  
  eta <- base_eta + m
  mu  <- 1/(1+exp(-eta))
  H   <- -sum(mu*(1-mu)) - invs2
  if(abs(H) < H_floor) H <- -H_floor
  v <- -1/H
  
  list(mean = m, var = v)
}



# ---- Data conversion helpers ----
# Convert M x n x n x T adjacency arrays to the edge-only M x E x T format
# used internally by the mixed-effects estimator.
adjacency_array_to_edges <- function(A, self_loop = FALSE) {
  if (length(dim(A)) != 4L) {
    stop("A must be an M x n x n x T array.", call. = FALSE)
  }
  M <- dim(A)[1]
  n <- dim(A)[2]
  if (dim(A)[3] != n) stop("The two node dimensions of A must be equal.", call. = FALSE)
  nT <- dim(A)[4]
  eg <- prepare_edge_index(n, self_loop = self_loop)
  Y <- array(0L, dim = c(M, eg$E, nT))
  for (t in seq_len(nT)) {
    for (m in seq_len(M)) {
      Y[m, , t] <- A[m, , , t][eg$idx]
    }
  }
  list(Y = Y, edge_index = eg)
}

# Convert edge-only data back to M x n x n x T symmetric adjacency arrays.
edges_to_adjacency_array <- function(Y, n, self_loop = FALSE) {
  if (length(dim(Y)) != 3L) stop("Y must be an M x E x T array.", call. = FALSE)
  M <- dim(Y)[1]
  nT <- dim(Y)[3]
  eg <- prepare_edge_index(n, self_loop = self_loop)
  if (dim(Y)[2] != eg$E) stop("The number of edge columns in Y is incompatible with n.", call. = FALSE)
  A <- array(0, dim = c(M, n, n, nT))
  for (t in seq_len(nT)) {
    for (m in seq_len(M)) {
      B <- matrix(0, n, n)
      B[eg$idx] <- Y[m, , t]
      B <- B + t(B)
      if (self_loop) diag(B) <- diag(B) / 2 else diag(B) <- 0
      A[m, , , t] <- B
    }
  }
  A
}

# ---- Output helpers shared by estimation / inference ----
parameter_metadata <- function(c) {
  c <- as.integer(factor(c, levels = sort(unique(c))))
  n <- length(c)
  G <- length(unique(c))

  theta_meta <- data.frame(
    parameter = paste0("theta_", seq_len(n)),
    type = "theta",
    node = seq_len(n),
    g = NA_integer_,
    l = NA_integer_,
    stringsAsFactors = FALSE
  )

  z_meta <- do.call(rbind, lapply(seq_len(G), function(l) {
    data.frame(
      parameter = paste0("Z_", seq_len(G), "_", l),
      type = "Z",
      node = NA_integer_,
      g = seq_len(G),
      l = l,
      stringsAsFactors = FALSE
    )
  }))

  sigma_meta <- data.frame(
    parameter = "sigma2",
    type = "sigma2",
    node = NA_integer_,
    g = NA_integer_,
    l = NA_integer_,
    stringsAsFactors = FALSE
  )

  rbind(theta_meta, z_meta, sigma_meta)
}

tidy_parameter_estimates <- function(fit, c, unique_Z = TRUE) {
  meta <- parameter_metadata(c)
  nT <- ncol(fit$theta)
  out <- vector("list", nT)

  for (t in seq_len(nT)) {
    est <- c(
      as.numeric(fit$theta[, t]),
      as.numeric(fit$Z[, , t]),
      fit$sigma[t]^2
    )
    d <- cbind(
      data.frame(time = t, u = t / nT),
      meta,
      estimate = est
    )
    if (unique_Z) {
      keep <- d$type != "Z" | d$g <= d$l
      d <- d[keep, , drop = FALSE]
    }
    rownames(d) <- NULL
    out[[t]] <- d
  }

  do.call(rbind, out)
}


# ---- Parallel worker helper ----
resolve_n_cores <- function(n_cores = NULL, n_tasks = 1L) {
  available <- suppressWarnings(parallel::detectCores(logical = TRUE))

  if (length(available) != 1L || is.na(available) ||
      !is.finite(available) || available < 1L) {
    available <- 1L
  }

  available <- as.integer(available)
  n_tasks <- max(1L, as.integer(n_tasks))

  if (is.null(n_cores)) {
    requested <- max(1L, available - 1L)
  } else {
    if (!is.numeric(n_cores) || length(n_cores) != 1L ||
        !is.finite(n_cores) || n_cores < 1) {
      stop("n_cores must be NULL or a positive integer.", call. = FALSE)
    }
    requested <- as.integer(n_cores)
  }

  min(requested, n_tasks, available)
}
