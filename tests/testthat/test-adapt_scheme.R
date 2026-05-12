library(graDiEnt)

obj <- function(x) sum((x - c(1, -1, 2))^2)
true_opt <- c(1, -1, 2)

cp_base <- function(scheme) {
  GetAlgoParams(
    n_params    = 3,
    n_iter      = 600,
    n_particles = 12,
    n_diff      = 2,
    init_sd     = 1,
    init_center = 0,
    adapt_scheme = scheme
  )
}

# ── each scheme converges ─────────────────────────────────────────────────────

test_that("adapt_scheme 'rand' converges to known optimum", {
  set.seed(42)
  out <- suppressMessages(optim_SQGDE(obj, cp_base("rand")))

  expect_true(all(is.finite(out$solution)))
  expect_lt(max(abs(out$solution - true_opt)), 0.1)
})

test_that("adapt_scheme 'best' converges to known optimum", {
  set.seed(42)
  out <- suppressMessages(optim_SQGDE(obj, cp_base("best")))

  expect_true(all(is.finite(out$solution)))
  expect_lt(max(abs(out$solution - true_opt)), 0.1)
})

test_that("adapt_scheme 'current' converges to known optimum", {
  set.seed(42)
  out <- suppressMessages(optim_SQGDE(obj, cp_base("current")))

  expect_true(all(is.finite(out$solution)))
  expect_lt(max(abs(out$solution - true_opt)), 0.1)
})

# ── each scheme returns correct output structure ──────────────────────────────

test_that("all schemes return valid list structure", {
  for (scheme in c("rand", "best", "current")) {
    set.seed(1)
    out <- suppressMessages(
      optim_SQGDE(obj, GetAlgoParams(n_params = 3, n_iter = 5,
                                     n_particles = 12, n_diff = 2,
                                     adapt_scheme = scheme))
    )
    expect_type(out, "list")
    expect_true(all(is.finite(out$solution)))
    expect_true(is.finite(out$weight))
    expect_length(out$solution, 3)
  }
})

# ── invalid scheme errors ─────────────────────────────────────────────────────

test_that("GetAlgoParams rejects invalid adapt_scheme", {
  expect_error(
    GetAlgoParams(n_params = 2, adapt_scheme = "gradient"),
    "invalid adaption scheme"
  )
  expect_error(
    GetAlgoParams(n_params = 2, adapt_scheme = "random"),
    "invalid adaption scheme"
  )
  expect_error(
    GetAlgoParams(n_params = 2, adapt_scheme = ""),
    "invalid adaption scheme"
  )
})
