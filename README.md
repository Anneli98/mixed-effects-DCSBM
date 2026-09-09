# Time-varying Mixed-Effects Degree-Corrected Stochastic Block Model

R code for estimating shared community structure and time-varying network parameters from replicated dynamic networks under a mixed-effects degree-corrected stochastic block model.

For an observed dataset, the main functions are:

1. `estimate_membership()` — estimate shared community memberships by residualized joint spectral embedding (RJSE), when memberships are unknown;

2. `select_bandwidth()` — select a common bandwidth by leave-one-out predictive marginal likelihood;

3. `mixed_effects_dcsbm()` — fit the time-varying mixed-effects DCSBM at a supplied bandwidth;

4. `mixed_effects_confint()` — construct pointwise confidence intervals for the fitted parameters.

All implementation details are stored in `R/subfunctions/`.

## Repository structure

```text
mixed-effects-DCSBM/
├── R/
│   ├── estimate_membership.R
│   ├── bandwidth_selection.R
│   ├── mixed_effects_dcsbm.R
│   ├── confidence_intervals.R
│   └── subfunctions/
│       ├── load_subfunctions.R
│       ├── base_functions.R
│       ├── laem_core.R
│       ├── marginal_likelihood.R
│       ├── inference_core.R
│       └── community/
│           ├── auxiliary_model.R
│           └── rjse_spectral_clustering.R
├── examples/
│   └── estimated_membership.R
├── README.md
├── .gitignore
├── .here
└── MixedEffectsDCSBM.Rproj
```

## Requirements

The code uses `here` for project-relative paths. The `parallel` package is distributed with R.

```r
install.packages(c("here", "statmod", "lme4"))
```

Community estimation by RJSE additionally requires the `networkSoSD` package.

For reliable path handling, open `MixedEffectsDCSBM.Rproj` before sourcing the main functions.

## Data format

The internal data format is an `M x E x T` array `Y`, where:

* `M` is the number of subjects;

* `T` is the number of time points;

* `E` is the number of undirected edges for a network.

Thus, `Y[m, , t]` contains the upper-triangular edge observations for subject `m` at time `t`.

If the data are stored as a full `M x n x n x T` adjacency array `A`, convert them by

```r
edge_data <- adjacency_array_to_edges(A)
Y <- edge_data$Y
```

## Quick start

```r
library(here)

source(here::here("R", "estimate_membership.R"))
source(here::here("R", "bandwidth_selection.R"))
source(here::here("R", "mixed_effects_dcsbm.R"))
source(here::here("R", "confidence_intervals.R"))
```

### 1. Community estimation

If the shared community memberships are unknown, estimate them by RJSE:

```r
community <- estimate_membership(
  Y = Y,
  G = 3,
  n_cores = 10,
  seed = 1
)

c_hat <- community$membership
```

If the memberships are already known, set `c_hat` directly to the observed label vector.

### 2. Bandwidth selection

Choose a candidate grid and select the bandwidth by leave-one-time-point-out cross-validation:

```r
bw <- select_bandwidth(
  Y = Y,
  c = c_hat,
  h_grid = h_grid,
  kernel = "epanechnikov",
  n_cores = 10
)

bw$h_opt
```

### 3. Fit the mixed-effects DCSBM

```r
fit <- mixed_effects_dcsbm(
  Y = Y,
  c = c_hat,
  h = bw$h_opt,
  kernel = "epanechnikov",
  n_cores = 10,
  max_iter = 1000,
  Q = 10,
  epsilon = 1e-5
)
```

The main fitted quantities are

```r
fit$theta       # n x T node effects
fit$Z           # G x G x T community-connectivity matrices
fit$sigma       # length-T random-effect standard deviations
fit$sigma2      # length-T random-effect variances
```

### 4. Confidence intervals

```r
ci <- mixed_effects_confint(
  Y = Y,
  c = c_hat,
  fit = fit,
  level = 0.95
)

ci$table
```

The confidence intervals use the tangent-space projected sandwich covariance, with optional smoothing-bias correction.

## Model

For subject `m`, node pair `(i,j)`, and time point `u^t`, the model is

```text
logit P{A_m,ij(u^t) = 1 | omega_m(u^t)}
  = theta_i(u^t)
  + theta_j(u^t)
  + Z_{c_i,c_j}(u^t)
  + omega_m(u^t),
```

with

```text
omega_m(u^t) ~ N(0, sigma^2(u^t)).
```

The community memberships are shared across time. The node effects satisfy the identifiability constraint

```text
sum_{i: c_i = g} theta_i(u^t) = 0,
```

and `Z(u^t)` is symmetric.

## Estimation

The parameter trajectories are estimated by a local approximate EM algorithm.

The main computational steps are:

* Laplace approximation for the subject-specific random intercept;

* Gaussian quadrature in the E-step;

* closed-form update for `sigma^2`;

* KKT-constrained one-step update for `theta`;

* one-step update and symmetrization for `Z`.

