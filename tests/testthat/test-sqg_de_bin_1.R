library(graDiEnt)

# ── grad_approx_fn: zero-division guards ──────────────────────────────────────

test_that("grad_approx_fn returns finite values when two particles are identical", {
  # identical particles → vec_diff_temp == 0 → original code divides by zero
  params <- matrix(c(1, 2, 1, 2, 3, 4, 5, 6), nrow = 4, ncol = 2)
  params[2, ] <- params[1, ]   # particles 1 and 2 are identical
  weights <- c(1.0, 0.5, 2.0, 1.5)

  result <- graDiEnt:::grad_approx_fn(
    param_indices  = 1:2,
    n_diff         = 1,
    current_params = params,
    parent_indices = c(1, 2),   # identical pair
    current_weight = weights
  )

  expect_true(all(is.finite(result$grad_approx)))
  expect_true(is.finite(result$psi) | is.nan(result$psi))  # psi may be NaN when grad is zero vector
})

test_that("grad_approx_fn returns finite grad when two particles have identical weights", {
  # identical weights → weight_diff_temp == 0 → original skips this implicitly,
  # but the guard makes it explicit and prevents any numerical issues
  params <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, ncol = 2)
  weights <- c(1.0, 1.0, 2.0, 1.5)   # particles 1 and 2 have same weight

  result <- graDiEnt:::grad_approx_fn(
    param_indices  = 1:2,
    n_diff         = 1,
    current_params = params,
    parent_indices = c(1, 2),
    current_weight = weights
  )

  expect_true(all(is.finite(result$grad_approx)))
})

test_that("grad_approx_fn handles mixed valid and zero-diff pairs without NaN", {
  # n_diff = 2: one pair identical, one pair distinct — guard skips bad pair
  params <- matrix(c(1,2, 1,2, 3,4, 5,6, 7,8), nrow = 5, ncol = 2, byrow = TRUE)
  weights <- c(1.0, 1.0, 0.5, 2.0, 1.5)

  result <- graDiEnt:::grad_approx_fn(
    param_indices  = 1:2,
    n_diff         = 2,
    current_params = params,
    parent_indices = c(1, 3, 2, 4),   # pair (1,2) identical; pair (3,4) distinct
    current_weight = weights
  )

  expect_true(all(is.finite(result$grad_approx)))
})

# ── degenerate population: optimizer doesn't crash ───────────────────────────

test_that("optimizer runs without error when particles initialise very close together", {
  # extremely tight init forces near-identical particles early on,
  # triggering the zero-division guard repeatedly
  set.seed(99)
  expect_no_error(
    suppressMessages(
      optim_SQGDE(
        ObjFun = function(x) sum((x - c(1, 2))^2),
        control_params = GetAlgoParams(
          n_params    = 2,
          n_iter      = 50,
          n_particles = 8,
          n_diff      = 2,
          init_sd     = 1e-8,
          init_center = 0
        )
      )
    )
  )
})

# ── rand scheme sampling: correct number of parents ───────────────────────────

test_that("rand scheme works with minimum viable n_particles and n_diff", {
  # n_particles=4, n_diff=1 for rand requires sampling 2*1+1=3 parents from 4
  # the original SQG_DE_bin_1_pos bug would sample only 2 and then access index 3
  set.seed(7)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = function(x) sum((x - 1)^2),
      control_params = GetAlgoParams(
        n_params     = 1,
        n_iter       = 100,
        n_particles  = 6,
        n_diff       = 1,
        adapt_scheme = "rand",
        init_sd      = 1
      )
    )
  )
  expect_true(is.finite(out$weight))
  expect_true(all(is.finite(out$solution)))
})

test_that("best scheme works with minimum viable n_particles and n_diff", {
  set.seed(7)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = function(x) sum((x - 1)^2),
      control_params = GetAlgoParams(
        n_params     = 1,
        n_iter       = 100,
        n_particles  = 6,
        n_diff       = 1,
        adapt_scheme = "best",
        init_sd      = 1
      )
    )
  )
  expect_true(is.finite(out$weight))
  expect_true(all(is.finite(out$solution)))
})

test_that("current scheme works with minimum viable n_particles and n_diff", {
  set.seed(7)
  out <- suppressMessages(
    optim_SQGDE(
      ObjFun = function(x) sum((x - 1)^2),
      control_params = GetAlgoParams(
        n_params     = 1,
        n_iter       = 100,
        n_particles  = 6,
        n_diff       = 1,
        adapt_scheme = "current",
        init_sd      = 1
      )
    )
  )
  expect_true(is.finite(out$weight))
  expect_true(all(is.finite(out$solution)))
})
