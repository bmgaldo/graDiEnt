library(graDiEnt)

# ── GetAlgoParams: bounds validation ─────────────────────────────────────────

test_that("GetAlgoParams rejects lower >= upper", {
  expect_error(GetAlgoParams(n_params = 2, lower = 1,      upper = 0),
               "lower must be strictly less than upper")
  expect_error(GetAlgoParams(n_params = 2, lower = c(0, 5), upper = c(1, 5)),
               "lower must be strictly less than upper")
})

test_that("GetAlgoParams rejects wrong-length bounds", {
  expect_error(GetAlgoParams(n_params = 2, lower = c(0, 0, 0)),
               "lower vector length must be 1 or n_params")
  expect_error(GetAlgoParams(n_params = 2, upper = c(1, 1, 1)),
               "upper vector length must be 1 or n_params")
})

test_that("GetAlgoParams stores bounds as length-n_params vectors", {
  cp <- GetAlgoParams(n_params = 3, lower = -1, upper = 2)
  expect_length(cp$lower, 3)
  expect_length(cp$upper, 3)
  expect_equal(cp$lower, c(-1, -1, -1))
  expect_equal(cp$upper, c(2,  2,  2))
})

test_that("GetAlgoParams accepts per-parameter bounds", {
  cp <- GetAlgoParams(n_params = 3, lower = c(-1, 0, -5), upper = c(1, 3, 5))
  expect_equal(cp$lower, c(-1,  0, -5))
  expect_equal(cp$upper, c( 1,  3,  5))
})

# ── optim_SQGDE: solution stays within symmetric scalar bounds ────────────────

test_that("solution stays within scalar bounds when optimum is inside", {
  set.seed(42)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = function(x) sum((x - c(1, 2))^2),
      control_params = GetAlgoParams(
        n_params    = 2,
        n_iter      = 400,
        n_particles = 10,
        n_diff      = 2,
        init_sd     = 1,
        init_center = 0,
        lower       = -5,
        upper       =  5
      )
    )
  )

  expect_true(all(out$solution >= -5))
  expect_true(all(out$solution <=  5))
  expect_lt(max(abs(out$solution - c(1, 2))), 0.15)
})

# ── optim_SQGDE: per-parameter asymmetric bounds ─────────────────────────────

test_that("solution stays within per-parameter asymmetric bounds", {
  set.seed(7)
  true_opt <- c(0.5, 2.5, -0.5)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = function(x) sum((x - true_opt)^2),
      control_params = GetAlgoParams(
        n_params    = 3,
        n_iter      = 500,
        n_particles = 12,
        n_diff      = 2,
        init_sd     = 0.3,
        init_center = c(0.5, 2, -0.5),
        lower       = c( 0,  1, -1),
        upper       = c( 1,  3,  0)
      )
    )
  )

  lower <- c(0, 1, -1)
  upper <- c(1, 3,  0)
  expect_true(all(out$solution >= lower))
  expect_true(all(out$solution <= upper))
  expect_lt(max(abs(out$solution - true_opt)), 0.15)
})

# ── optim_SQGDE: optimum outside bounds → clips to boundary ──────────────────

test_that("solution clips to boundary when true optimum is outside bounds", {
  set.seed(99)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = function(x) sum((x - c(-10, -10))^2),
      control_params = GetAlgoParams(
        n_params    = 2,
        n_iter      = 400,
        n_particles = 10,
        n_diff      = 2,
        init_sd     = 0.5,
        init_center = 1,
        lower       = 0,
        upper       = 2
      )
    )
  )

  expect_true(all(out$solution >= 0))
  expect_true(all(out$solution <= 2))
  # constrained optimum is at lower boundary (0, 0)
  expect_lt(max(abs(out$solution - c(0, 0))), 0.1)
})

# ── optim_SQGDE: unbounded by default (no behaviour change) ──────────────────

test_that("default bounds (-Inf / Inf) do not alter unbounded behaviour", {
  obj <- function(x) sum((x - c(-2, 3))^2)

  set.seed(5)
  out_nobound <- suppressMessages(
    optim_SQGDE(
      ObjFun = obj,
      control_params = GetAlgoParams(
        n_params    = 2,
        n_iter      = 400,
        n_particles = 10,
        n_diff      = 2,
        init_sd     = 1
      )
    )
  )

  set.seed(5)
  out_explicitinf <- suppressMessages(
    optim_SQGDE(
      ObjFun = obj,
      control_params = GetAlgoParams(
        n_params    = 2,
        n_iter      = 400,
        n_particles = 10,
        n_diff      = 2,
        init_sd     = 1,
        lower       = -Inf,
        upper       =  Inf
      )
    )
  )

  expect_equal(out_nobound$solution, out_explicitinf$solution)
})
