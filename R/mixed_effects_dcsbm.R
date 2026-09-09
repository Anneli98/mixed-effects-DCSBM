# ============================================================
# mixed_effects_dcsbm.R
# User-facing estimator for the time-varying mixed-effects DCSBM
# ============================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. Install it with install.packages('here').",
    call. = FALSE
  )
}

here::i_am("R/mixed_effects_dcsbm.R")

source(
  here::here(
    "R",
    "subfunctions",
    "load_subfunctions.R"
  )
)


#' Fit the time-varying mixed-effects DCSBM
#'
#' @param Y An M x E x T edge-only array.
#' @param c Length-n vector of known community memberships.
#' @param h Positive bandwidth for temporal smoothing.
#' @param kernel Kernel name: "epanechnikov", "uniform",
#'   "triangular", or "gaussian". A kernel function may also be supplied.
#' @param n_cores Number of workers. If NULL, choose automatically.
#' @param self_loop Whether self-loops are included.
#' @param max_iter Maximum number of LA-EM iterations.
#' @param Q Number of Gaussian quadrature points.
#' @param epsilon Convergence tolerance.
#' @param damp Damping factor for the fixed-effect updates.
#' @param verbose Whether to print iteration information.
#' @param omega_newton_maxit Maximum Newton iterations for omega.
#' @param omega_newton_tol Newton tolerance for omega.
#'
#' @return A fitted mixed-effects DCSBM object.
mixed_effects_dcsbm <- function(
    Y,
    c,
    h,
    kernel = "epanechnikov",
    n_cores = NULL,
    self_loop = FALSE,
    max_iter = 1000,
    Q = 10,
    epsilon = 1e-5,
    damp = 0.8,
    verbose = FALSE,
    omega_newton_maxit = 20,
    omega_newton_tol = 1e-5
) {

  .medcsbm_require("statmod")
  .medcsbm_require("lme4")

  if (length(dim(Y)) != 3L) {
    stop(
      "Y must be an M x E x T edge-only array.",
      call. = FALSE
    )
  }

  if (
    anyNA(Y) ||
      !all(Y == 0 | Y == 1)
  ) {
    stop(
      "Y must contain only binary edge observations 0/1.",
      call. = FALSE
    )
  }

  if (
    anyNA(c) ||
      length(c) < 2L
  ) {
    stop(
      "c must contain one non-missing community label for each node.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(h) ||
      length(h) != 1L ||
      !is.finite(h) ||
      h <= 0
  ) {
    stop(
      "h must be a positive finite number.",
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
  nT <- dim(Y)[3]
  G <- length(unique(c))

  eg <- prepare_edge_index(
    n,
    self_loop = self_loop
  )

  if (eg$E != dim(Y)[2]) {
    stop(
      sprintf(
        "Y has E=%d edge columns, but n=%d implies E=%d.",
        dim(Y)[2],
        n,
        eg$E
      ),
      call. = FALSE
    )
  }

  ker <- get_kernel(kernel)

  n_workers <- resolve_n_cores(
    n_cores,
    n_tasks = nT
  )

  # Internal numerical starting values.
  .start_values <- function(Y_t) {

    M <- nrow(Y_t)

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

  theta0 <- matrix(
    0,
    nrow = n,
    ncol = nT
  )

  Z0 <- array(
    NA_real_,
    dim = c(G, G, nT)
  )

  sigma0 <- numeric(nT)

  for (tt in seq_len(nT)) {

    init_t <- .start_values(
      Y[
        ,
        ,
        tt,
        drop = FALSE
      ][, , 1]
    )

    theta0[, tt] <-
      center_theta_by_community(
        init_t$theta,
        c
      )

    Z0[, , tt] <-
      symmetrize_Z(
        init_t$Z
      )

    sigma0[tt] <-
      max(
        as.numeric(init_t$sigma),
        1e-8
      )
  }

  fit <- fit_laem_path(
    Y = Y,
    c = c,
    h = h,
    kernel = ker$fun,
    self_loop = self_loop,
    smooth = TRUE,
    leave_one_time_out = FALSE,
    max_iter = max_iter,
    Q = Q,
    epsilon = epsilon,
    damp = damp,
    verbose = verbose,

    # Q is always monitored and used for the original
    # convergence safeguard / best-Q return rule.
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

  fit$settings$kernel <- ker$name
  fit$settings$kernel_support <- ker$support

  fit$point_estimates <-
    tidy_parameter_estimates(
      fit,
      c,
      unique_Z = TRUE
    )

  class(fit) <- unique(
    c(
      "mixed_effects_dcsbm_fit",
      class(fit)
    )
  )

  fit
}


print.mixed_effects_dcsbm_fit <- function(x, ...) {

  cat(
    "Time-varying mixed-effects DCSBM fit\n"
  )

  cat(
    sprintf(
      "  h: %g\n",
      x$settings$h
    )
  )

  cat(
    sprintf(
      "  kernel: %s\n",
      x$settings$kernel
    )
  )

  cat(
    sprintf(
      "  workers: %d\n",
      x$settings$n_cores
    )
  )

  cat(
    sprintf(
      "  theta: %d x %d\n",
      nrow(x$theta),
      ncol(x$theta)
    )
  )

  cat(
    sprintf(
      "  Z: %d x %d x %d\n",
      dim(x$Z)[1],
      dim(x$Z)[2],
      dim(x$Z)[3]
    )
  )

  cat(
    sprintf(
      "  sigma: length %d\n",
      length(x$sigma)
    )
  )

  invisible(x)
}
