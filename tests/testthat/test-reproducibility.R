library(graDiEnt)

# ── Reproducibility: identical seeds in the same session ─────────────────────

test_that("same-seed runs return identical results (including traces)", {
  set.seed(42)
  n  <- 200
  x  <- rnorm(n)
  y  <- 0.5 + 1.5 * x + rnorm(n, sd = 0.5)

  ols_obj <- function(params, x, y) {
    sum((y - params[1] - params[2] * x)^2)
  }

  cp <- GetAlgoParams(
    n_params     = 2,
    n_iter       = 100,
    n_particles  = 12,
    n_diff       = 2,
    init_center  = 0,
    init_sd     = 1,
    return_trace = TRUE
  )

  set.seed(123)
  out1 <- suppressMessages(
    optim_SQGDE(ObjFun = ols_obj, control_params = cp, x = x, y = y)
  )

  set.seed(123)
  out2 <- suppressMessages(
    optim_SQGDE(ObjFun = ols_obj, control_params = cp, x = x, y = y)
  )

  expect_identical(out1, out2)
})

test_that("different seeds produce different results", {
  set.seed(42)
  n  <- 200
  x  <- rnorm(n)
  y  <- 0.5 + 1.5 * x + rnorm(n, sd = 0.5)

  ols_obj <- function(params, x, y) {
    sum((y - params[1] - params[2] * x)^2)
  }

  cp <- GetAlgoParams(
    n_params    = 2,
    n_iter      = 100,
    n_particles = 12,
    n_diff      = 2,
    init_center = 0,
    init_sd     = 1
  )

  set.seed(123)
  out1 <- suppressMessages(
    optim_SQGDE(ObjFun = ols_obj, control_params = cp, x = x, y = y)
  )

  set.seed(456)
  out2 <- suppressMessages(
    optim_SQGDE(ObjFun = ols_obj, control_params = cp, x = x, y = y)
  )

  expect_false(identical(out1, out2))
})