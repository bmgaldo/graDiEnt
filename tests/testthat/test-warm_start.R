library(graDiEnt)

obj_2d <- function(x) sum((x - c(1, 2))^2)

# ── output always contains last_particles and last_weights ────────────────────

test_that("optim_SQGDE returns last_particles and last_weights without return_trace", {
  set.seed(1)
  out <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 8, n_diff = 2))
  )

  expect_true("last_particles" %in% names(out))
  expect_true("last_weights"   %in% names(out))
  expect_equal(dim(out$last_particles), c(8L, 2L))
  expect_length(out$last_weights, 8)
  expect_true(all(is.finite(out$last_particles)))
  expect_true(all(is.finite(out$last_weights)))
})

test_that("optim_SQGDE returns last_particles and last_weights with return_trace", {
  set.seed(1)
  out <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 8,
                              n_diff = 2, return_trace = TRUE))
  )

  expect_true("last_particles" %in% names(out))
  expect_true("last_weights"   %in% names(out))
  expect_equal(dim(out$last_particles), c(8L, 2L))
  expect_length(out$last_weights, 8)
})

# ── warm start produces valid output ─────────────────────────────────────────

test_that("warm start returns valid solution with correct structure", {
  set.seed(1)
  out1 <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 100, n_particles = 8,
                              n_diff = 2, init_sd = 1))
  )

  set.seed(2)
  out2 <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 100, n_particles = 8, n_diff = 2),
                warm_start = out1)
  )

  expect_type(out2, "list")
  expect_named(out2, c("solution", "weight", "last_particles", "last_weights", "converged"))
  expect_length(out2$solution, 2)
  expect_true(all(is.finite(out2$solution)))
  expect_true(is.finite(out2$weight))
})

# ── warm start does not degrade solution ─────────────────────────────────────

test_that("warm start weight is no worse than cold start weight", {
  set.seed(42)
  out1 <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 200, n_particles = 10,
                              n_diff = 2, init_sd = 1))
  )

  out2 <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 200, n_particles = 10, n_diff = 2),
                warm_start = out1)
  )

  expect_lte(out2$weight, out1$weight)
})

# ── chained warm starts converge to solution ─────────────────────────────────

test_that("chained warm starts converge closer to true optimum", {
  set.seed(7)
  cp <- GetAlgoParams(n_params = 2, n_iter = 150, n_particles = 10,
                      n_diff = 2, init_sd = 2)

  out <- suppressMessages(optim_SQGDE(obj_2d, cp))
  for (i in 1:3) {
    out <- suppressMessages(
      optim_SQGDE(obj_2d,
                  GetAlgoParams(n_params = 2, n_iter = 150, n_particles = 10, n_diff = 2),
                  warm_start = out)
    )
  }

  expect_lt(max(abs(out$solution - c(1, 2))), 0.1)
})

# ── warm start works with bounds ──────────────────────────────────────────────

test_that("warm start respects bounds from new control_params", {
  set.seed(5)
  out1 <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 100, n_particles = 8,
                              n_diff = 2, lower = -5, upper = 5))
  )

  out2 <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 100, n_particles = 8,
                              n_diff = 2, lower = -5, upper = 5),
                warm_start = out1)
  )

  expect_true(all(out2$solution >= -5))
  expect_true(all(out2$solution <=  5))
  expect_lte(out2$weight, out1$weight)
})

# ── validation errors ─────────────────────────────────────────────────────────

test_that("warm_start missing last_particles or last_weights throws error", {
  bad <- list(solution = c(0, 0), weight = 1)
  expect_error(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 8, n_diff = 2),
                warm_start = bad),
    "missing last_particles or last_weights"
  )
})

test_that("warm_start n_particles mismatch throws error", {
  set.seed(1)
  out1 <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 8, n_diff = 2))
  )

  expect_error(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 12, n_diff = 2),
                warm_start = out1),
    "control_params expects"
  )
})

test_that("warm_start n_params mismatch throws error", {
  set.seed(1)
  out1 <- suppressMessages(
    optim_SQGDE(obj_2d,
                GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 8, n_diff = 2))
  )

  expect_error(
    optim_SQGDE(function(x) sum((x - 1:3)^2),
                GetAlgoParams(n_params = 3, n_iter = 50, n_particles = 8, n_diff = 2),
                warm_start = out1),
    "control_params expects"
  )
})
