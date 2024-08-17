library(testthat)

test_that("GetAlgoParams assigns defaults and handles input correctly", {

  # Test with minimal input
  params <- GetAlgoParams(n_params = 10)

  expect_equal(params$n_params, 10)
  expect_equal(params$n_particles, 30)  # 3 * n_params
  expect_equal(params$n_iter, 1000)
  expect_equal(params$init_sd, 0.01)
  expect_equal(params$init_center, 0)
  expect_equal(params$n_cores_use, 1)
  expect_equal(params$step_size, 2.38 / sqrt(2 * 10))
  expect_equal(params$jitter_size, 1e-6)
  expect_equal(params$crossover_rate, 1)
  expect_equal(params$parallel_type, 'none')
  expect_equal(params$return_trace, FALSE)
  expect_equal(params$thin, 1)
  expect_equal(params$purify, Inf)
  expect_equal(params$adapt_scheme, 'rand')
  expect_equal(params$give_up_init, 100)
  expect_equal(params$stop_tol, 1e-4)
  expect_equal(params$stop_check, 10)
  expect_equal(params$converge_crit, 'stdev')
  expect_equal(params$var_list, NULL)
  expect_equal(params$message_int, 100)

  # Test with custom input
  params <- GetAlgoParams(
    n_params = 5,
    n_particles = 20,
    n_iter = 500,
    init_sd = 0.05,
    init_center = 1,
    n_cores_use = 4,
    step_size = 1.5,
    jitter_size = 1e-5,
    crossover_rate = 0.8,
    parallel_type = 'PSOCK',
    return_trace = TRUE,
    thin = 5,
    purify = 10,
    adapt_scheme = 'best',
    give_up_init = 50,
    stop_tol = 1e-5,
    stop_check = 5,
    converge_crit = 'percent',
    var_list = c("var1", "var2"),
    message_int = 50
  )

  expect_equal(params$n_params, 5)
  expect_equal(params$n_particles, 20)
  expect_equal(params$n_iter, 500)
  expect_equal(params$init_sd, 0.05)
  expect_equal(params$init_center, 1)
  expect_equal(params$n_cores_use, 4)
  expect_equal(params$step_size, 1.5)
  expect_equal(params$jitter_size, 1e-5)
  expect_equal(params$crossover_rate, 0.8)
  expect_equal(params$parallel_type, 'PSOCK')
  expect_equal(params$return_trace, TRUE)
  expect_equal(params$thin, 5)
  expect_equal(params$purify, 10)
  expect_equal(params$adapt_scheme, 'best')
  expect_equal(params$give_up_init, 50)
  expect_equal(params$stop_tol, 1e-5)
  expect_equal(params$stop_check, 5)
  expect_equal(params$converge_crit, 'percent')
  expect_equal(params$var_list, c("var1", "var2"))
  expect_equal(params$message_int, 50)

  # Test for error handling
  expect_error(GetAlgoParams(n_params = -1), "n_params must be a postitive integer scalar")
  expect_error(GetAlgoParams(n_params = 10, n_particles = 2), "n_particles must be a postitive integer scalar, and atleast 4")
  expect_error(GetAlgoParams(n_params = 10, step_size = -1), "step_size must be positive and real-valued")
  expect_error(GetAlgoParams(n_params = 10, crossover_rate = 1.5), "crossover_rate must be a numeric scalar on the interval \\(0,1\\]")
  expect_error(GetAlgoParams(n_params = 10, parallel_type = "invalid"), "invalid parallel_type")
  expect_error(GetAlgoParams(n_params = 10, adapt_scheme = "invalid"), "invalid adaption scheme")
  expect_error(GetAlgoParams(n_params = 10, thin = 0), "thin must be a scalar postive integer")
  expect_error(GetAlgoParams(n_params = 10, purify = -1), "purify must be a positive integer or Inf")
})

