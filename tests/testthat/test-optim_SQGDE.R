library(graDiEnt)

# ── MLE: factorizable multivariate normal (known sigma, estimate means) ──────

test_that("MLE factorizable MVN recovers sample means", {
  set.seed(42)
  true_mu <- c(-1, 1, 0, 2)
  data <- matrix(rnorm(500 * 4, mean = true_mu, sd = 1), nrow = 500, ncol = 4)
  analytic_mu <- colMeans(data)

  neg_log_lik <- function(x, data) {
    -sum(
      dnorm(data[, 1], x[1], sd = 1, log = TRUE) +
      dnorm(data[, 2], x[2], sd = 1, log = TRUE) +
      dnorm(data[, 3], x[3], sd = 1, log = TRUE) +
      dnorm(data[, 4], x[4], sd = 1, log = TRUE)
    )
  }

  set.seed(123)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = neg_log_lik,
      control_params = GetAlgoParams(
        n_params    = 4,
        n_iter      = 1000,
        n_particles = 12,
        n_diff      = 2,
        init_center = 0,
        init_sd     = 2
      ),
      data = data
    )
  )

  expect_type(out, "list")
  expect_named(out, c("solution", "weight", "converged"))
  expect_length(out$solution, 4)
  expect_true(all(is.finite(out$solution)))
  expect_lt(max(abs(out$solution - analytic_mu)), 0.1)
})

# ── Univariate objective (n_params = 1) ─────────────────────────────────────

test_that("univariate objective recovers sample mean", {
  set.seed(42)
  obs <- rnorm(300, mean = 2.5, sd = 1)
  analytic_mu <- mean(obs)

  neg_log_lik_1d <- function(x, obs) {
    -sum(dnorm(obs, mean = x[1], sd = 1, log = TRUE))
  }

  set.seed(123)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = neg_log_lik_1d,
      control_params = GetAlgoParams(
        n_params    = 1,
        n_iter      = 500,
        n_particles = 6,
        n_diff      = 1,
        init_center = 0,
        init_sd     = 2
      ),
      obs = obs
    )
  )

  expect_type(out, "list")
  expect_named(out, c("solution", "weight", "converged"))
  expect_length(out$solution, 1)
  expect_true(is.finite(out$solution))
  expect_lt(abs(out$solution - analytic_mu), 0.1)
})

# ── OLS: univariate regression (intercept + 1 slope) ────────────────────────

test_that("OLS univariate regression recovers lm coefficients", {
  set.seed(42)
  n  <- 200
  x  <- rnorm(n)
  y  <- 0.5 + 1.5 * x + rnorm(n, sd = 0.5)

  lm_coef <- unname(coef(lm(y ~ x)))

  ols_obj <- function(params, x, y) {
    sum((y - params[1] - params[2] * x)^2)
  }

  set.seed(123)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = ols_obj,
      control_params = GetAlgoParams(
        n_params    = 2,
        n_iter      = 1000,
        n_particles = 12,
        n_diff      = 2,
        init_center = 0,
        init_sd     = 1
      ),
      x = x,
      y = y
    )
  )

  expect_type(out, "list")
  expect_length(out$solution, 2)
  expect_true(all(is.finite(out$solution)))
  expect_lt(max(abs(out$solution - lm_coef)), 0.1)
})

# ── OLS: multivariate regression (intercept + 3 slopes) ─────────────────────

test_that("OLS multivariate regression recovers lm coefficients", {
  set.seed(42)
  n         <- 300
  Xmat      <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  true_beta <- c(0.3, -0.5, 1.0, 0.7)
  y         <- cbind(1, Xmat) %*% true_beta + rnorm(n, sd = 0.5)

  lm_coef <- unname(coef(lm(y ~ Xmat)))

  ols_obj <- function(params, Xmat, y) {
    sum((y - params[1] - Xmat %*% params[2:4])^2)
  }

  set.seed(123)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = ols_obj,
      control_params = GetAlgoParams(
        n_params    = 4,
        n_iter      = 1000,
        n_particles = 12,
        n_diff      = 2,
        init_center = 0,
        init_sd     = 1
      ),
      Xmat = Xmat,
      y    = y
    )
  )

  expect_type(out, "list")
  expect_length(out$solution, 4)
  expect_true(all(is.finite(out$solution)))
  expect_lt(max(abs(out$solution - lm_coef)), 0.1)
})

# ── Parallel: PSOCK (2 cores) ────────────────────────────────────────────────

test_that("PSOCK parallel execution returns valid solution", {
  skip_on_cran()

  set.seed(42)
  true_mu <- c(-1, 1, 0, 2)
  data_par <- matrix(rnorm(500 * 4, mean = true_mu, sd = 1), nrow = 500, ncol = 4)
  analytic_mu <- colMeans(data_par)

  neg_log_lik <- function(x, data_par) {
    -sum(
      dnorm(data_par[, 1], x[1], sd = 1, log = TRUE) +
      dnorm(data_par[, 2], x[2], sd = 1, log = TRUE) +
      dnorm(data_par[, 3], x[3], sd = 1, log = TRUE) +
      dnorm(data_par[, 4], x[4], sd = 1, log = TRUE)
    )
  }

  set.seed(123)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = neg_log_lik,
      control_params = GetAlgoParams(
        n_params      = 4,
        n_iter        = 500,
        n_particles   = 12,
        n_diff        = 2,
        init_center   = 0,
        init_sd       = 2,
        parallel_type = "PSOCK",
        n_cores_use   = 2
      ),
      data_par = data_par
    )
  )

  expect_type(out, "list")
  expect_named(out, c("solution", "weight", "converged"))
  expect_length(out$solution, 4)
  expect_true(all(is.finite(out$solution)))
  expect_lt(max(abs(out$solution - analytic_mu)), 0.1)
})

# ── Parallel: FORK (2 cores, non-Windows only) ───────────────────────────────

test_that("FORK parallel execution returns valid solution", {
  skip_on_cran()
  skip_on_os("windows")

  set.seed(42)
  true_mu <- c(-1, 1, 0, 2)
  data_par <- matrix(rnorm(500 * 4, mean = true_mu, sd = 1), nrow = 500, ncol = 4)
  analytic_mu <- colMeans(data_par)

  neg_log_lik <- function(x, data_par) {
    -sum(
      dnorm(data_par[, 1], x[1], sd = 1, log = TRUE) +
      dnorm(data_par[, 2], x[2], sd = 1, log = TRUE) +
      dnorm(data_par[, 3], x[3], sd = 1, log = TRUE) +
      dnorm(data_par[, 4], x[4], sd = 1, log = TRUE)
    )
  }

  set.seed(123)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = neg_log_lik,
      control_params = GetAlgoParams(
        n_params      = 4,
        n_iter        = 1000,
        n_particles   = 12,
        n_diff        = 2,
        init_center   = 0,
        init_sd       = 2,
        parallel_type = "FORK",
        n_cores_use   = 2
      ),
      data_par = data_par
    )
  )

  expect_type(out, "list")
  expect_named(out, c("solution", "weight", "converged"))
  expect_length(out$solution, 4)
  expect_true(all(is.finite(out$solution)))
  expect_lt(max(abs(out$solution - analytic_mu)), 0.1)
})
