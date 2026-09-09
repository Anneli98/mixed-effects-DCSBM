# ============================================================
# Pointwise confidence intervals for the mixed-effects DCSBM
# ============================================================
# Uses the tangent-space projected sandwich covariance from the
# asymptotic theory. Optional O(h^2) bias correction is included.

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}

here::i_am("R/confidence_intervals.R")
source(here::here("R", "subfunctions", "load_subfunctions.R"))

mixed_effects_confint <- function(Y, c, fit,
                         times = NULL,
                         level = 0.95,
                         kernel = NULL,
                         bias_correct = TRUE,
                         Q = 15,
                         ridge = 1e-8,
                         omega_newton_maxit = 50,
                         omega_newton_tol = 1e-8,
                         center_scores = TRUE,
                         unique_Z = TRUE,
                         keep_details = FALSE) {
  if (is.null(times)) times <- seq_len(dim(Y)[3])
  times <- sort(unique(as.integer(times)))
  times <- times[times >= 1L & times <= dim(Y)[3]]
  if (length(times) == 0L) stop("No valid time indices supplied.", call. = FALSE)

  if (isTRUE(fit$settings$leave_one_time_out)) {
    warning(
      "The supplied fit is a leave-one-time-out fit. For final inference, use the full-sample fit ",
      "with leave_one_time_out = FALSE.",
      call. = FALSE
    )
  }

  if (is.null(kernel)) kernel <- fit$settings$kernel
  ker <- get_kernel(kernel)
  support <- fit$settings$kernel_support
  if (is.null(support)) support <- ker$support

  meta <- parameter_metadata(c)
  nT <- dim(Y)[3]
  tables <- vector("list", length(times))
  details <- if (keep_details) vector("list", length(times)) else NULL

  for (kk in seq_along(times)) {
    t <- times[kk]

    v <- estimate_variance_time(
      Y = Y, c = c, t = t, fit = fit,
      h = fit$settings$h,
      kernel = ker$fun,
      self_loop = isTRUE(fit$settings$self_loop),
      smooth = TRUE,
      Q = Q,
      kernel_support = support,
      omega_newton_maxit = omega_newton_maxit,
      omega_newton_tol = omega_newton_tol,
      ridge = ridge,
      center_scores = center_scores,
      return_vcov = keep_details
    )

    b <- if (bias_correct) {
      estimate_smoothing_bias(v, fit, t = t, h = fit$settings$h,
               kernel = ker$fun, support = support)$bias
    } else {
      rep(0, length(v$se_xi))
    }

    ci <- construct_confidence_interval(fit, t, bias = b, se = v$se_xi, level = level)
    d <- cbind(
      data.frame(time = t, u = t / nT),
      meta,
      as.data.frame(ci)
    )

    if (unique_Z) {
      keep <- d$type != "Z" | d$g <= d$l
      d <- d[keep, , drop = FALSE]
    }
    rownames(d) <- NULL
    tables[[kk]] <- d
    if (keep_details) details[[kk]] <- v
  }

  out <- list(
    table = do.call(rbind, tables),
    details = details,
    level = level,
    bias_correct = bias_correct,
    h = fit$settings$h,
    kernel = ker$name,
    times = times
  )
  class(out) <- c("mixed_effects_confint", class(out))
  out
}
