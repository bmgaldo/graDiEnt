test_that("ObjFun_list validation", {
  # Valid input: single function
  valid_single_function <- function(x) {
    x^2
  }

  expect_message(optim_SQGDE(ObjFun_list = valid_single_function,
                             control_params = GetAlgoParams(n_params = 1)))

  # Valid input: list of functions
  valid_list_functions <- list(
    function(x) {x^2},
    function(x) {x^3}
  )
  expect_message(optim_SQGDE(ObjFun_list = valid_list_functions,
                             control_params =
                               GetAlgoParams(n_params = 1,
                                             param_ind_to_update_list =
                                               list(1, 1))))

  # Invalid input: non-function
  invalid_non_function <- "not a function"
  expect_error(optim_SQGDE(ObjFun_list = invalid_non_function,
                           control_params = GetAlgoParams(n_params = 1)))

  # Invalid input: list with non-function
  invalid_list_with_non_function <- list(function(x) {x^2}, "not a function")
  expect_error(optim_SQGDE(ObjFun_list = invalid_list_with_non_function,
                           control_params = GetAlgoParams(n_params = 1)))

  # Invalid input: inconsistent lengths
  inconsistent_length_list <- list(function(x) {x^2}, function(x) {x^3})
  expect_error(optim_SQGDE(ObjFun_list = inconsistent_length_list,
                           control_params = GetAlgoParams(n_params = 1)))
})

test_that("Parallel cluster initialization", {
  # Mock control_params
  parallel_type = "PSOCK"
  n_cores_use = 2
  parallel_seed = 123
  x <<- 1
  y <<- 2
  varlist = c("x", "y")

  # Test PSOCK cluster initialization
  # debugonce(optim_SQGDE)
  expect_message(optim_SQGDE(ObjFun_list = function(x) {x^2},
                            control_params =
                              GetAlgoParams(n_params = 1,
                                parallel_type = parallel_type,
                                n_cores_use = n_cores_use,
                                parallel_seed = parallel_seed,
                                varlist = varlist)
                            )
                )

  # Test FORK cluster initialization (if supported on the system)
  if(Sys.info()['sysname'] != "Windows"){
    parallel_type <- "FORK"
    expect_message(optim_SQGDE(ObjFun_list = function(x) {x^2},
                                control_params =
                                  GetAlgoParams(n_params = 1,
                                                parallel_type = parallel_type,
                                                n_cores_use = n_cores_use,
                                                parallel_seed = parallel_seed,
                                                varlist = varlist)
    ))
  }
  # Test seed setting
  # A more indirect approach could be to run the optimization multiple times with the same seed and
  # check if the results are consistent.
  set.seed(43210)
  output_1 <- optim_SQGDE(ObjFun_list = function(x) {x^2},
                             control_params =
                               GetAlgoParams(n_params = 1,
                                             parallel_type = parallel_type,
                                             n_cores_use = n_cores_use,
                                             parallel_seed = parallel_seed,
                                             varlist = varlist)
  )

  set.seed(43210)
  output_2 <- optim_SQGDE(ObjFun_list = function(x) {x^2},
                          control_params =
                            GetAlgoParams(n_params = 1,
                                          parallel_type = parallel_type,
                                          n_cores_use = n_cores_use,
                                          parallel_seed = parallel_seed,
                                          varlist = varlist)
  )
  expect_equal(output_1$solution, output_2$solution)
})


test_that("Percent convergence", {
  # Mock control_params
  parallel_type = "PSOCK"
  n_cores_use = 2
  parallel_seed = 123
  converge_crit <- "percent"

  # Test PSOCK cluster initialization
  set.seed(123)
  expect_message(optim_SQGDE(ObjFun_list = function(x) {x^2},
                             control_params =
                               GetAlgoParams(n_params = 1,
                                             parallel_type = parallel_type,
                                             n_cores_use = n_cores_use,
                                             parallel_seed = parallel_seed,
                                             converge_crit = converge_crit
                                              )
  )
  )
})

test_that("stdev convergence", {
  # Mock control_params
  parallel_type = "PSOCK"
  n_cores_use = 2
  parallel_seed = 123
  converge_crit <- "stdev"

  # Test PSOCK cluster initialization
  set.seed(1234)
  expect_message(optim_SQGDE(ObjFun_list = function(x) {x^2},
                             control_params =
                               GetAlgoParams(n_params = 1,
                                             parallel_type = parallel_type,
                                             n_cores_use = n_cores_use,
                                             parallel_seed = parallel_seed,
                                             converge_crit = converge_crit,
                                             stop_tol = 1e-1
                               )
  )
  )
})
