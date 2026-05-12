library(graDiEnt)

ALL_NAMES <- c("n_params", "n_particles", "n_iter", "init_sd", "init_center",
               "lower", "upper", "bounds_type", "n_cores_use", "step_size",
               "crossover_rate", "jitter_size", "parallel_type", "recovery_path",
               "recovery_freq", "thin", "purify", "n_iters_per_particle",
               "return_trace", "n_diff", "adapt_scheme", "give_up_init",
               "stop_tol", "stop_check", "converge_crit", "trace_print_freq")

# ── output structure ──────────────────────────────────────────────────────────

test_that("GetAlgoParams returns list with all expected names", {
  cp <- GetAlgoParams(n_params = 3)
  expect_named(cp, ALL_NAMES, ignore.order = TRUE)
})

# ── default values ────────────────────────────────────────────────────────────

test_that("n_particles default is max(3*n_params, 4)", {
  expect_equal(GetAlgoParams(n_params = 5)$n_particles, 15L)
  expect_equal(GetAlgoParams(n_params = 1)$n_particles, 4L)  # max(3,4)=4
})

test_that("n_iter default is 1000", {
  expect_equal(GetAlgoParams(n_params = 2)$n_iter, 1000L)
})

test_that("step_size default is 2.38/sqrt(2*n_params)", {
  cp <- GetAlgoParams(n_params = 4)
  expect_equal(cp$step_size, 2.38 / sqrt(2 * 4))
})

test_that("adapt_scheme default is 'rand'", {
  expect_equal(GetAlgoParams(n_params = 2)$adapt_scheme, "rand")
})

test_that("bounds_type default is 'reflect'", {
  expect_equal(GetAlgoParams(n_params = 2)$bounds_type, "reflect")
})

test_that("lower/upper default to vectors of -Inf/Inf length n_params", {
  cp <- GetAlgoParams(n_params = 3)
  expect_equal(cp$lower, rep(-Inf, 3))
  expect_equal(cp$upper, rep(Inf, 3))
})

test_that("n_iters_per_particle = floor(n_iter / thin)", {
  cp <- GetAlgoParams(n_params = 2, n_iter = 100, thin = 3)
  expect_equal(cp$n_iters_per_particle, floor(100 / 3))
})

test_that("converge_crit default is 'stdev'", {
  expect_equal(GetAlgoParams(n_params = 2)$converge_crit, "stdev")
})

test_that("jitter_size default is 1e-6", {
  expect_equal(GetAlgoParams(n_params = 2)$jitter_size, 1e-6)
})

test_that("crossover_rate default is 1", {
  expect_equal(GetAlgoParams(n_params = 2)$crossover_rate, 1)
})

test_that("purify default is Inf", {
  expect_equal(GetAlgoParams(n_params = 2)$purify, Inf)
})

test_that("recovery_path default is NULL", {
  expect_null(GetAlgoParams(n_params = 2)$recovery_path)
})

test_that("recovery_freq default is 1", {
  expect_equal(GetAlgoParams(n_params = 2)$recovery_freq, 1L)
})

# ── n_params validation ───────────────────────────────────────────────────────

test_that("n_params is coerced to integer", {
  cp <- GetAlgoParams(n_params = 3.0)
  expect_identical(cp$n_params, 3L)
})

test_that("n_params < 1 errors", {
  expect_error(GetAlgoParams(n_params = 0), "n_params")
})

test_that("n_params non-finite errors", {
  expect_error(GetAlgoParams(n_params = Inf), "n_params")
})

test_that("n_params length > 1 errors", {
  expect_error(GetAlgoParams(n_params = c(2, 3)), "scalar")
})

# ── n_particles validation ────────────────────────────────────────────────────

test_that("n_particles < 4 errors", {
  expect_error(GetAlgoParams(n_params = 2, n_particles = 3), "n_particles")
})

test_that("n_particles non-finite errors", {
  expect_error(GetAlgoParams(n_params = 2, n_particles = Inf), "n_particles")
})

test_that("n_particles coerced to integer", {
  cp <- GetAlgoParams(n_params = 2, n_particles = 8.0)
  expect_identical(cp$n_particles, 8L)
})

# ── n_iter validation ─────────────────────────────────────────────────────────

test_that("n_iter < 4 errors", {
  expect_error(GetAlgoParams(n_params = 2, n_iter = 3), "n_iter")
})

test_that("n_iter coerced to integer", {
  expect_identical(GetAlgoParams(n_params = 2, n_iter = 100.0)$n_iter, 100L)
})

# ── n_diff validation ─────────────────────────────────────────────────────────

test_that("n_diff > n_particles/2 errors", {
  expect_error(
    GetAlgoParams(n_params = 2, n_particles = 6, n_diff = 4),
    "n_diff"
  )
})

test_that("n_diff < 1 errors", {
  expect_error(GetAlgoParams(n_params = 2, n_diff = 0), "n_diff")
})

test_that("n_diff coerced to integer", {
  expect_identical(GetAlgoParams(n_params = 2, n_diff = 2.0)$n_diff, 2L)
})

# ── init_sd validation ────────────────────────────────────────────────────────

test_that("init_sd <= 0 errors", {
  expect_error(GetAlgoParams(n_params = 2, init_sd = 0), "init_sd")
  expect_error(GetAlgoParams(n_params = 2, init_sd = -1), "init_sd")
})

test_that("init_sd non-finite errors", {
  expect_error(GetAlgoParams(n_params = 2, init_sd = Inf), "init_sd")
})

test_that("init_sd wrong length errors", {
  expect_error(GetAlgoParams(n_params = 2, init_sd = c(0.1, 0.2, 0.3)), "init_sd")
})

test_that("init_sd length n_params accepted", {
  expect_no_error(GetAlgoParams(n_params = 3, init_sd = c(0.1, 0.2, 0.3)))
})

# ── init_center validation ────────────────────────────────────────────────────

test_that("init_center non-finite errors", {
  expect_error(GetAlgoParams(n_params = 2, init_center = Inf), "init_center")
})

test_that("init_center wrong length errors", {
  expect_error(GetAlgoParams(n_params = 2, init_center = c(0, 0, 0)), "init_center")
})

test_that("init_center length n_params accepted", {
  expect_no_error(GetAlgoParams(n_params = 3, init_center = c(1, 2, 3)))
})

# ── step_size validation ──────────────────────────────────────────────────────

test_that("step_size <= 0 errors", {
  expect_error(GetAlgoParams(n_params = 2, step_size = 0), "step_size")
  expect_error(GetAlgoParams(n_params = 2, step_size = -1), "step_size")
})

test_that("step_size non-finite errors", {
  expect_error(GetAlgoParams(n_params = 2, step_size = Inf), "step_size")
})

test_that("step_size length > 1 errors", {
  expect_error(GetAlgoParams(n_params = 2, step_size = c(0.5, 0.5)), "step_size")
})

# ── jitter_size validation ────────────────────────────────────────────────────

test_that("jitter_size = 0 accepted (turns off jitter)", {
  expect_no_error(GetAlgoParams(n_params = 2, jitter_size = 0))
})

test_that("jitter_size < 0 errors", {
  expect_error(GetAlgoParams(n_params = 2, jitter_size = -1e-6), "jitter_size")
})

test_that("jitter_size non-finite errors", {
  expect_error(GetAlgoParams(n_params = 2, jitter_size = Inf), "jitter_size")
})

test_that("jitter_size length > 1 errors", {
  expect_error(GetAlgoParams(n_params = 2, jitter_size = c(1e-6, 1e-6)), "jitter_size")
})

# ── crossover_rate validation ─────────────────────────────────────────────────

test_that("crossover_rate = 0 accepted", {
  expect_no_error(GetAlgoParams(n_params = 2, crossover_rate = 0))
})

test_that("crossover_rate = 1 accepted", {
  expect_no_error(GetAlgoParams(n_params = 2, crossover_rate = 1))
})

test_that("crossover_rate > 1 errors", {
  expect_error(GetAlgoParams(n_params = 2, crossover_rate = 1.1), "crossover_rate")
})

test_that("crossover_rate < 0 errors", {
  expect_error(GetAlgoParams(n_params = 2, crossover_rate = -0.1), "crossover_rate")
})

# ── n_cores_use validation ────────────────────────────────────────────────────

test_that("n_cores_use < 1 errors", {
  expect_error(GetAlgoParams(n_params = 2, n_cores_use = 0), "n_cores_use")
})

test_that("n_cores_use coerced to integer", {
  expect_identical(GetAlgoParams(n_params = 2, n_cores_use = 2.0)$n_cores_use, 2L)
})

# ── parallel_type validation ──────────────────────────────────────────────────

test_that("invalid parallel_type errors", {
  expect_error(GetAlgoParams(n_params = 2, parallel_type = "MPI"), "parallel_type")
})

test_that("valid parallel_type values accepted", {
  expect_no_error(GetAlgoParams(n_params = 2, parallel_type = "none"))
  expect_no_error(GetAlgoParams(n_params = 2, parallel_type = "PSOCK"))
  expect_no_error(GetAlgoParams(n_params = 2, parallel_type = "FORK"))
})

# ── thin validation ───────────────────────────────────────────────────────────

test_that("thin < 1 errors", {
  expect_error(GetAlgoParams(n_params = 2, thin = 0), "thin")
})

test_that("thin coerced to integer", {
  expect_identical(GetAlgoParams(n_params = 2, thin = 5.0)$thin, 5L)
})

# ── purify validation ─────────────────────────────────────────────────────────

test_that("purify = Inf accepted", {
  expect_no_error(GetAlgoParams(n_params = 2, purify = Inf))
})

test_that("purify positive integer accepted", {
  expect_no_error(GetAlgoParams(n_params = 2, purify = 10))
})

test_that("purify < 1 errors", {
  expect_error(GetAlgoParams(n_params = 2, purify = 0), "purify")
})

# ── give_up_init validation ───────────────────────────────────────────────────

test_that("give_up_init < 1 errors", {
  expect_error(GetAlgoParams(n_params = 2, give_up_init = 0), "give_up_init")
})

# ── stop_check validation ─────────────────────────────────────────────────────

test_that("stop_check < 2 errors", {
  expect_error(GetAlgoParams(n_params = 2, stop_check = 1), "stop_check")
})

# ── stop_tol validation ───────────────────────────────────────────────────────

test_that("stop_tol < 0 errors", {
  expect_error(GetAlgoParams(n_params = 2, stop_tol = -1e-5), "stop_tol")
})

test_that("stop_tol = 0 accepted", {
  expect_no_error(GetAlgoParams(n_params = 2, stop_tol = 0))
})

# ── converge_crit validation ──────────────────────────────────────────────────

test_that("invalid converge_crit errors", {
  expect_error(GetAlgoParams(n_params = 2, converge_crit = "max"), "converge_crit")
})

test_that("valid converge_crit values accepted", {
  expect_no_error(GetAlgoParams(n_params = 2, converge_crit = "stdev"))
  expect_no_error(GetAlgoParams(n_params = 2, converge_crit = "percent"))
})
