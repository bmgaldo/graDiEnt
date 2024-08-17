library(testthat)

# Define a simple objective function for testing
test_obj_fun <- function(x) {
  return(sum((x - 3)^2))
}

test_that("optim_SQGDE runs without error and returns expected structure", {

  # Basic test with minimal parameters
  control_params <- GetAlgoParams(n_params = 2, n_iter = 10, n_particles = 6)
  result <- optim_SQGDE(ObjFun = test_obj_fun, control_params = control_params)

  # Check that the result has the expected structure
  expect_type(result, "list")
  expect_true(all(c("solution", "weight", "converged") %in% names(result)))
  expect_length(result$solution, control_params$n_params)
  expect_type(result$weight, "double")
  expect_type(result$converged, "logical")

  # Check that the solution is numeric and of the correct length
  expect_type(result$solution, "double")
  expect_equal(length(result$solution), control_params$n_params)
})

test_that("optim_SQGDE can handle more iterations and particles", {

  # Test with more iterations and particles
  control_params <- GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 20)
  result <- optim_SQGDE(ObjFun = test_obj_fun, control_params = control_params)

  # Check that the solution is closer to the known minimum
  expect_true(all(result$solution > 2 & result$solution < 4))
})

test_that("optim_SQGDE handles convergence correctly", {

  # Define a control parameter setup with a convergence criterion
  control_params <- GetAlgoParams(n_params = 2, n_iter = 100, n_particles = 10, stop_tol = 1e-6, stop_check = 10, converge_crit = "stdev")
  result <- optim_SQGDE(ObjFun = test_obj_fun, control_params = control_params)

  # Check if the function converges
  expect_true(result$converged)
})

test_that("optim_SQGDE handles parallelization correctly", {

  # Test with parallelization (PSOCK)
  control_params <- GetAlgoParams(n_params = 2, n_iter = 20, n_particles = 10, parallel_type = "PSOCK", n_cores_use = 2)
  result <- optim_SQGDE(ObjFun = test_obj_fun, control_params = control_params)

  # Check that the result has the expected structure
  expect_type(result, "list")
  expect_true(all(c("solution", "weight", "converged") %in% names(result)))
  expect_length(result$solution, control_params$n_params)
})

test_that("optim_SQGDE returns particle traces when requested", {

  # Test with return_trace = TRUE
  control_params <- GetAlgoParams(n_params = 2, n_iter = 30, n_particles = 10, return_trace = TRUE)
  result <- optim_SQGDE(ObjFun = test_obj_fun, control_params = control_params)

  # Check that the particle traces are returned
  expect_type(result$particles_trace, "double")
  expect_type(result$weights_trace, "double")
  expect_equal(dim(result$particles_trace), c(control_params$n_iters_per_particle, control_params$n_particles, control_params$n_params))
  expect_equal(dim(result$weights_trace), c(control_params$n_iters_per_particle, control_params$n_particles))
})

test_that("optim_SQGDE throws error on initialization failure", {

  # Test with an objective function that always returns Inf
  bad_obj_fun <- function(x) Inf
  control_params <- GetAlgoParams(n_params = 2, n_iter = 10, n_particles = 5, give_up_init = 3)

  expect_error(optim_SQGDE(ObjFun = bad_obj_fun, control_params = control_params), "population initialization failed")
})
