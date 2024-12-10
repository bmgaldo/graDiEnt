
n_params <- 1


# Test basic functionality with valid arguments
test_that("GetAlgoParams works with valid arguments", {
  params <- GetAlgoParams(n_params = 5)
  expect_equal(length(params), 27)  # Check all arguments are returned
  expect_equal(params$n_params, 5)
  expect_equal(params$n_particles, 15)  # Default value
})

# Test cases for n_params validation
test_that("n_params validation works as expected", {
  # Valid input
  expect_silent(GetAlgoParams(n_params = 5))

  # Invalid input types
  expect_error(GetAlgoParams(n_params = "5"), "ERROR: n_params is not a finite number")
  expect_error(GetAlgoParams(n_params = TRUE), "ERROR: n_params is not a finite number")
  expect_error(GetAlgoParams(n_params = c(1, 2)), "ERROR: n_params is not a finite number")

  # Invalid numeric values
  expect_error(GetAlgoParams(n_params = -1), "ERROR: n_params must be a postitive integer scalar")
  expect_error(GetAlgoParams(n_params = 0), "ERROR: n_params must be a postitive integer scalar")
  expect_error(GetAlgoParams(n_params = Inf), "ERROR: n_params is not a finite number")
  expect_error(GetAlgoParams(n_params = NaN), "ERROR: n_params is not a finite number")
})

test_that("param_ind_to_update_list validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, param_ind_to_update_list = NULL))

  # Valid input: list of logical vectors
  expect_silent(GetAlgoParams(n_params = 5, param_ind_to_update_list = list(c(TRUE, FALSE, TRUE, FALSE, TRUE))))

  # Valid input: list of numeric vectors
  expect_error(GetAlgoParams(n_params = 5,
                             param_ind_to_update_list = list(c(1, 3, 5))))

  # Invalid input: non-list
  expect_error(GetAlgoParams(n_params = 5,
                             param_ind_to_update_list = c(1, 2, 3)))

  # Invalid input: list with wrong length
  expect_error(GetAlgoParams(n_params = 5,
                             param_ind_to_update_list =
                               list(c(TRUE, FALSE), c(TRUE, FALSE))))

  # Invalid input: list with non-logical or non-numeric elements
  expect_error(GetAlgoParams(n_params = 5,
                             param_ind_to_update_list = list("a", 2)))
  expect_error(GetAlgoParams(n_params = 5,
                             param_ind_to_update_list = list(c(1, 2, 3, 4))))
})


test_that("resample_weight validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, resample_weight = NULL))

  # Valid input: TRUE
  expect_silent(GetAlgoParams(n_params = 5, resample_weight = TRUE))

  # Valid input: FALSE
  expect_silent(GetAlgoParams(n_params = 5, resample_weight = FALSE))

  # Invalid input: non-logical
  expect_silent(GetAlgoParams(n_params = 5, resample_weight = "TRUE"))
  expect_silent(GetAlgoParams(n_params = 5, resample_weight = 1))
})

test_that("n_particles validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, n_particles = NULL))

  # Valid input: integer greater than or equal to 4
  expect_silent(GetAlgoParams(n_params = 5, n_particles = 10))

  # Invalid input: non-integer
  expect_error(GetAlgoParams(n_params = 5, n_particles = 2.5),
               "ERROR: n_particles must be a postitive integer scalar")

  # Invalid input: integer less than 4
  expect_error(
    GetAlgoParams(n_params = 5, n_particles = 3),
    "ERROR: n_particles must be a postitive integer scalar, and atleast 4")

  # Invalid input: non-finite
  expect_error(GetAlgoParams(n_params = 5, n_particles = Inf),
               "ERROR: n_particles is not finite") |>
    expect_warning()
  expect_error(GetAlgoParams(n_params = 5, n_particles = NaN),
               "ERROR: n_particles is not finite")
})

test_that("n_iter validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, n_iter = NULL))

  # Valid input: integer greater than or equal to 4
  expect_silent(GetAlgoParams(n_params = 5, n_iter = 100))

  # Invalid input: non-integer
  expect_error(GetAlgoParams(n_params = 5, n_iter = 2.5),
               "ERROR: n_iter must be a postitive integer scalar")

  # Invalid input: integer less than 4
  expect_error(GetAlgoParams(n_params = 5, n_iter = 3),
               "ERROR: n_iter must be a postitive integer scalar, and atleast 4")

  # Invalid input: non-finite
  expect_error(GetAlgoParams(n_params = 5, n_iter = Inf),
               "ERROR: n_iter is not finite") |>
    expect_warning()
  expect_error(GetAlgoParams(n_params = 5, n_iter = NaN),
               "ERROR: n_iter is not finite")
})

test_that("init_sd validation works as expected", {
  # Valid input: single positive number
  expect_silent(GetAlgoParams(n_params = 5, init_sd = 0.1))

  # Valid input: vector of positive numbers with length equal to n_params
  expect_silent(GetAlgoParams(n_params = 5,
                              init_sd = c(0.1, 0.2, 0.3, 0.4, 0.5)))

  # Invalid input: non-numeric
  expect_silent(GetAlgoParams(n_params = 5, init_sd = "0.1"))

  # Invalid input: negative value
  expect_error(GetAlgoParams(n_params = 5, init_sd = -0.1),
               "ERROR: init_sd must be positive and real-valued")

  # Invalid input: complex number
  expect_error(GetAlgoParams(n_params = 5, init_sd = 1i),
               "ERROR: init_sd must be positive and real-valued") |>
    expect_warning()

  # Invalid input: vector with incorrect length
  expect_error(GetAlgoParams(n_params = 5, init_sd = c(0.1, 0.2)),
               "ERROR: init_sd vector length must be 1 or n_params")

  # Invalid input: non-finite value
  expect_error(GetAlgoParams(n_params = 5, init_sd = Inf),
               "ERROR: init_sd is not finite")
  expect_error(GetAlgoParams(n_params = 5, init_sd = NaN),
               "ERROR: init_sd is not finite")
})

test_that("init_center validation works as expected", {
  # Valid input: single number
  expect_silent(GetAlgoParams(n_params = 5, init_center = 0))

  # Valid input: vector of numbers with length equal to n_params
  expect_silent(GetAlgoParams(n_params = 5, init_center = c(0, 1, 2, 3, 4)))

  # Invalid input: non-numeric
  expect_silent(GetAlgoParams(n_params = 5, init_center = "0"))

  # Invalid input: complex number
  expect_warning(GetAlgoParams(n_params = 5, init_center = 1i))

  # Invalid input: vector with incorrect length
  expect_error(GetAlgoParams(n_params = 5, init_center = c(0, 1)),
               "ERROR: init_center vector length must be 1 or n_params")

  # Invalid input: non-finite value
  expect_error(GetAlgoParams(n_params = 5, init_center = Inf),
               "ERROR: init_center is not finite")
  expect_error(GetAlgoParams(n_params = 5, init_center = NaN),
               "ERROR: init_center is not finite")
})

test_that("n_cores_use validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, n_cores_use = NULL))

  # Valid input: integer greater than or equal to 1
  expect_silent(GetAlgoParams(n_params = 5, n_cores_use = 4))

  # Invalid input: non-integer
  expect_silent(GetAlgoParams(n_params = 5, n_cores_use = 2.5))

  # Invalid input: integer less than 1
  expect_error(
    GetAlgoParams(n_params = 5, n_cores_use = 0),
    "ERROR: n_cores_use must be a postitive integer scalar, and atleast 1")

  # Invalid input: non-finite
  expect_error(GetAlgoParams(n_params = 5, n_cores_use = Inf),
               "ERROR: n_cores_use is not finite") |>
    expect_warning()
  expect_error(GetAlgoParams(n_params = 5, n_cores_use = NaN),
               "ERROR: n_cores_use is not finite")
})

test_that("step_size validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, step_size = NULL))

  # Valid input: positive number
  expect_silent(GetAlgoParams(n_params = 5, step_size = 0.5))

  # Invalid input: non-numeric
  expect_error(GetAlgoParams(n_params = 5, step_size = "0.5"),
               "ERROR: step_size is not finite")

  # Invalid input: negative number
  expect_error(GetAlgoParams(n_params = 5, step_size = -0.5),
               "ERROR: step_size must be positive and real-valued")

  # Invalid input: complex number
  expect_error(GetAlgoParams(n_params = 5, step_size = 1i))

  # Invalid input: vector
  expect_error(GetAlgoParams(n_params = 5, step_size = c(0.5, 0.6)),
               "ERROR: step_size vector length must be 1 ")

  # Invalid input: non-finite value
  expect_error(GetAlgoParams(n_params = 5, step_size = Inf),
               "ERROR: step_size is not finite")
  expect_error(GetAlgoParams(n_params = 5, step_size = NaN),
               "ERROR: step_size is not finite")
})

test_that("jitter_size validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, jitter_size = NULL))

  # Valid input: positive number
  expect_silent(GetAlgoParams(n_params = 5, jitter_size = 0.1))

  # Invalid input: non-numeric
  expect_error(GetAlgoParams(n_params = 5, jitter_size = "0.1"),
               "ERROR: jitter_size is not finite")

  # Invalid input: negative number
  expect_error(GetAlgoParams(n_params = 5, jitter_size = -0.1),
               "ERROR: jitter_size must be positive and real-valued")

  # Invalid input: complex number
  expect_error(GetAlgoParams(n_params = 5, jitter_size = 1i))

  # Invalid input: vector
  expect_error(GetAlgoParams(n_params = 5, jitter_size = c(0.1, 0.2)),
               "ERROR: jitter_size vector length must be 1 ")

  # Invalid input: non-finite value
  expect_error(GetAlgoParams(n_params = 5, jitter_size = Inf),
               "ERROR: jitter_size is not finite")
  expect_error(GetAlgoParams(n_params = 5, jitter_size = NaN),
               "ERROR: jitter_size is not finite")
})

test_that("crossover_rate validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, crossover_rate = NULL))

  # Valid input: number between 0 and 1
  expect_silent(GetAlgoParams(n_params = 5, crossover_rate = 0.5))

  # Invalid input: non-numeric
  expect_silent(GetAlgoParams(n_params = 5, crossover_rate = "0.5"))

  # Invalid input: negative number
  expect_error(GetAlgoParams(n_params = 5, crossover_rate = -0.5))

  # Invalid input: number greater than 1
  expect_error(GetAlgoParams(n_params = 5, crossover_rate = 1.5))

  # Invalid input: complex number
  expect_error(GetAlgoParams(n_params = 5, crossover_rate = 1i)) |>
    expect_warning()

  # Invalid input: vector
  expect_error(GetAlgoParams(n_params = 5, crossover_rate = c(0.5, 0.6)))

  # Invalid input: non-finite value
  expect_error(GetAlgoParams(n_params = 5, crossover_rate = Inf),
               "ERROR: crossover_rate is not finite")
  expect_error(GetAlgoParams(n_params = 5, crossover_rate = NaN),
               "ERROR: crossover_rate is not finite")
})

test_that("parallel_type validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, parallel_type = NULL))

  # Valid input: "none"
  expect_silent(GetAlgoParams(n_params = 5, parallel_type = "none"))

  # Valid input: "FORK"
  expect_silent(GetAlgoParams(n_params = 5, parallel_type = "FORK"))

  # Valid input: "PSOCK"
  expect_silent(GetAlgoParams(n_params = 5, parallel_type = "PSOCK"))

  # Invalid input: non-character
  expect_error(GetAlgoParams(n_params = 5, parallel_type = 1),
               "ERROR: invalid parallel_type.")

  # Invalid input: invalid string
  expect_error(GetAlgoParams(n_params = 5, parallel_type = "invalid"),
               "ERROR: invalid parallel_type.")
})

test_that("parallel_seed validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, parallel_seed = NULL))

  # Valid input: positive integer
  expect_silent(GetAlgoParams(n_params = 5, parallel_seed = 123))

  # Invalid input: non-integer
  expect_error(GetAlgoParams(n_params = 5, parallel_seed = 1.23),
               "ERROR: invalid parallel_seed.")

  # Invalid input: negative integer
  expect_error(GetAlgoParams(n_params = 5, parallel_seed = -1),
               "ERROR: invalid parallel_seed.")
})

test_that("converge_crit validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, converge_crit = NULL))

  # Valid input: "stdev"
  expect_silent(GetAlgoParams(n_params = 5, converge_crit = "stdev"))

  # Valid input: "percent"
  expect_silent(GetAlgoParams(n_params = 5, converge_crit = "percent"))

  # Invalid input: non-character
  expect_error(GetAlgoParams(n_params = 5, converge_crit = 1),
               "ERROR: invalid converge_crit.")

  # Invalid input: invalid string
  expect_error(GetAlgoParams(n_params = 5, converge_crit = "invalid"),
               "ERROR: invalid converge_crit.")
})

test_that("adapt_scheme validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, adapt_scheme = NULL))

  # Valid input: "rand"
  expect_silent(GetAlgoParams(n_params = 5, adapt_scheme = "rand"))

  # Valid input: "current"
  expect_silent(GetAlgoParams(n_params = 5, adapt_scheme = "current"))

  # Valid input: "best"
  expect_silent(GetAlgoParams(n_params = 5, adapt_scheme = "best"))

  # Invalid input: non-character
  expect_error(GetAlgoParams(n_params = 5, adapt_scheme = 1),
               "ERROR: invalid adaption scheme.")

  # Invalid input: invalid string
  expect_error(GetAlgoParams(n_params = 5, adapt_scheme = "invalid"),
               "ERROR: invalid adaption scheme.")
})

test_that("thin validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, thin = NULL))

  # Valid input: positive integer
  expect_silent(GetAlgoParams(n_params = 5, thin = 10))


  expect_silent(GetAlgoParams(n_params = 5, thin = 2.5))

  # Invalid input: negative integer
  expect_error(GetAlgoParams(n_params = 5, thin = -1),
               "ERROR: thin must be a scalar postive integer")

  # Invalid input: non-finite value
  expect_error(GetAlgoParams(n_params = 5, thin = Inf),
               "ERROR: thin is not finite") |>
    expect_warning()
  expect_error(GetAlgoParams(n_params = 5, thin = NaN),
               "ERROR: thin is not finite")
})

test_that("n_iters_per_particle calculation and validation works as expected", {
  # Valid input: n_iter > thin
  expect_silent(GetAlgoParams(n_params = 5, n_iter = 100, thin = 10))

  # Valid input: n_iter == thin
  expect_silent(GetAlgoParams(n_params = 5, n_iter = 10, thin = 10))

  # Invalid input: n_iter < thin
  expect_error(
    GetAlgoParams(n_params = 5, n_iter = 5, thin = 10),
    "ERROR: number of stored particle value is negative or non finite.")
})

test_that("purify validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, purify = NULL))

  # Valid input: positive integer
  expect_silent(GetAlgoParams(n_params = 5, purify = 10))

  # Valid input: Inf
  expect_warning(
    GetAlgoParams(n_params = 5, purify = Inf),
    'Warning: purify was not finite. Defaulting to not using Purify.')

  # Decimal input: non-integer coercion
  expect_silent(GetAlgoParams(n_params = 5, purify = 2.5))

  # Invalid input: negative integer
  expect_error(GetAlgoParams(n_params = 5, purify = -1),
               "ERROR: purify must be a positive integer or Inf")

  # Invalid input: non-finite value (other than Inf)
  expect_warning(
    GetAlgoParams(n_params = 5, purify = NaN),
    'Warning: purify was not finite. Defaulting to not using Purify.')
})

test_that("n_diff validation works as expected", {
  # Valid input: positive integer less than or equal to n_particles/2
  expect_silent(GetAlgoParams(n_params = 5, n_particles = 10, n_diff = 2))

  # Decimal input: non-integer coercion
  expect_silent(GetAlgoParams(n_params = 5, n_particles = 10, n_diff = 2.5))

  # Invalid input: negative integer
  expect_error(GetAlgoParams(n_params = 5, n_particles = 10, n_diff = -1),
               "ERROR: n_diff must be a scalar postive integer")

  # Invalid input: integer greater than n_particles/2
  expect_error(GetAlgoParams(n_params = 5, n_particles = 10, n_diff = 6),
               "ERROR: n_diff cannot exceed n_particles/2")

  # Invalid input: non-finite value
  expect_error(GetAlgoParams(n_params = 5, n_particles = 10, n_diff = Inf),
               "ERROR: n_diff is not finite") |>
    expect_warning()
  expect_error(GetAlgoParams(n_params = 5, n_particles = 10, n_diff = NaN),
               "ERROR: n_diff is not finite")
})

test_that("give_up_init validation works as expected", {
  # Valid input: NULL
  expect_warning(
    GetAlgoParams(n_params = 5, give_up_init = NULL),
    'Warning: give_up_init was null or not finite. Using default value of 100.')

  # Valid input: positive integer
  expect_silent(GetAlgoParams(n_params = 5, give_up_init = 10))

  # Decimal input: non-integer coercion
  expect_silent(GetAlgoParams(n_params = 5, give_up_init = 2.5))

  # Invalid input: negative integer
  expect_error(GetAlgoParams(n_params = 5, give_up_init = -1),
               "ERROR: give_up_init must be a scalar positive integer")

  # Invalid input: non-finite value
  expect_warning(
    GetAlgoParams(n_params = 5, give_up_init = Inf),
    'Warning: give_up_init was null or not finite. Using default value of 100.')
  expect_warning(
    GetAlgoParams(n_params = 5, give_up_init = NaN),
    'Warning: give_up_init was null or not finite. Using default value of 100.')
})

test_that("stop_check validation works as expected", {
  # Valid input: NULL
  expect_warning(
    GetAlgoParams(n_params = 5, stop_check = NULL),
    'Warning: stop_check was null or non-finite. Using default value of 10.')

  # Valid input: positive integer greater than 2
  expect_silent(GetAlgoParams(n_params = 5, stop_check = 10))

  # Invalid input: rounding decimal
  expect_silent(GetAlgoParams(n_params = 5, stop_check = 2.5))

  # Invalid input: integer less than or equal to 2
  expect_error(
    GetAlgoParams(n_params = 5, stop_check = 1),
    "ERROR: stop_check must be a scalar positive integer and greater than 1")
  expect_error(
    GetAlgoParams(n_params = 5, stop_check = 1),
    "ERROR: stop_check must be a scalar positive integer and greater than 1")
  expect_error(
    GetAlgoParams(n_params = 5, stop_check = 0),
    "ERROR: stop_check must be a scalar positive integer and greater than 1")
  expect_error(
    GetAlgoParams(n_params = 5, stop_check = -1),
    "ERROR: stop_check must be a scalar positive integer and greater than 1")

  # Invalid input: non-finite value
  expect_warning(
    GetAlgoParams(n_params = 5,
                  stop_check = Inf),
    'Warning: stop_check was null or non-finite. Using default value of 10.')
  expect_warning(
    GetAlgoParams(n_params = 5, stop_check = NaN),
    'Warning: stop_check was null or non-finite. Using default value of 10.')
})

test_that("stop_tol validation works as expected", {
  # Valid input: single positive scalar
  expect_silent(GetAlgoParams(n_params = 5, stop_tol = 1e-4))

  # character coercion: non-numeric
  expect_silent(GetAlgoParams(n_params = 5, stop_tol = "0.1"))

  # Invalid input: negative value
  expect_error(GetAlgoParams(n_params = 5, stop_tol = -0.1),
               "ERROR: stop_tol must be a scalar positive")

  # Invalid input: vector
  expect_error(GetAlgoParams(n_params = 5, stop_tol = c(1e-4, 2e-3)),
               "ERROR: stop_tol must be a scalar positive")

  # Invalid input: NULL
  # This should set a default value
  expect_silent(GetAlgoParams(n_params = 5, stop_tol = NULL))

  # Invalid input: NA
  # This should set a default value
  expect_silent(GetAlgoParams(n_params = 5, stop_tol = NA))
})

test_that("varlist validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, varlist = NULL))

  # Valid input: list of strings
  expect_silent(GetAlgoParams(n_params = 5, varlist = list("x", "y", "z")))

  # Edge input: single variable
  expect_silent(GetAlgoParams(n_params = 5, varlist = "x"))

  # Invalid input: list with non-string elements
  expect_error(
    GetAlgoParams(
      n_params = 5,
      varlist = list("x", 1, "z")),
    "ERROR: each element of stop varlist must be a string of variable and function name!")
})

test_that("save_int validation works as expected", {
  # Valid input: NULL
  expect_silent(GetAlgoParams(n_params = 5, save_int = NULL))

  # Valid input: positive integer
  expect_silent(GetAlgoParams(n_params = 5, save_int = 10))

  # Valid input: negative integer
  expect_silent(GetAlgoParams(n_params = 5, save_int = -1))

  # Invalid input: non-integer
  expect_silent(GetAlgoParams(n_params = 5,
                              save_int = 2.5))

  # Invalid input: non-finite value
  expect_error(GetAlgoParams(n_params = 5, save_int = Inf),
               "ERROR: save_int is not finite")
  expect_error(GetAlgoParams(n_params = 5, save_int = NaN),
               "ERROR: save_int is not finite")

  # Invalid input: zero
  expect_error(GetAlgoParams(n_params = 5,
                             save_int = 0))
})

test_that("save_rds_string validation works as expected", {
  # Valid input: character string ending in .rds
  expect_silent(GetAlgoParams(n_params = 5,
                              save_int = 10,
                              save_rds_string = "my_model.rds"))

  # Invalid input: non-character
  expect_error(
    GetAlgoParams(n_params = 5, save_int = 10, save_rds_string = 123),
    "ERROR: save_rds_string must be a single CHARACTER sting that ends in .rds")

  # Invalid input: character string not ending in .rds
  expect_warning(GetAlgoParams(n_params = 5,
                               save_int = 10,
                               save_rds_string = "my_model"))

  # Invalid input: vector of strings
  expect_error(
    GetAlgoParams(n_params = 5,
                  save_int = 10,
                  save_rds_string = c("model1.rds", "model2.rds")),
    "ERROR: save_rds_string must be a SINGLE character sting that ends in .rds")
})

test_that("outfile_string validation works as expected", {
  # Valid input: character string ending in .rds
  expect_silent(GetAlgoParams(n_params = 5,
                              outfile_string = "my_par_model.txt"))

  # Invalid input: non-character
  expect_error(
    GetAlgoParams(n_params = 5, outfile_string = 123),
    "ERROR: outfile_string must be a single CHARACTER sting that ends in .rds")

  # Invalid input: character string not ending in .rds
  expect_warning(GetAlgoParams(n_params = 5,
                               outfile_string = "my_par_model"))

  # Invalid input: vector of strings
  expect_error(
    GetAlgoParams(n_params = 5,
                  outfile_string = c("my_par_model1.txt", "my_par_model2.txt")),
    "ERROR: outfile_string must be a SINGLE character sting that ends in .rds")
})

# Test for errors with invalid input types
test_that("GetAlgoParams handles invalid input types", {
  expect_error(GetAlgoParams(n_params = "invalid"))  # Non-numeric n_params
  expect_error(GetAlgoParams(n_params,
                             n_particles = c(1, 2, 3)))  # Non-scalar n_particles
  expect_error(GetAlgoParams(n_params,
                             init_sd = -1))  # Negative init_sd
  # Crossover rate out of range
  expect_error(GetAlgoParams(n_params,
                             crossover_rate = 1.1))
})

# Test for invalid values in list arguments
test_that("GetAlgoParams handles invalid values in list arguments", {
  # Non-logical element
  expect_error(GetAlgoParams(n_params,
                             param_ind_to_update_list = list("FALSE")))
})

# Test for invalid parallel type
test_that("GetAlgoParams handles invalid parallel type", {
  # Not a valid option
  expect_error(GetAlgoParams(n_params,
                             parallel_type = "invalid_type"))
})

# Test for invalid purify value
test_that("GetAlgoParams handles invalid purify value", {
  expect_error(GetAlgoParams(n_params,
                             purify = 0))  # Needs to be positive or Inf
})

# Test for invalid n_diff value
test_that("GetAlgoParams handles invalid n_diff value", {
  expect_error(GetAlgoParams(n_params,
                             n_diff = 0))  # Needs to be positive
  expect_error(GetAlgoParams(n_params,
                             n_diff = 6))  # Cannot exceed n_particles/2
})

# Test for invalid give_up_init value
test_that("GetAlgoParams handles invalid give_up_init value", {
  expect_error(GetAlgoParams(n_params,
                             give_up_init = 0))  # Needs to be positive
})

# Test for invalid stop_check value
test_that("GetAlgoParams handles invalid stop_check value", {
  expect_error(GetAlgoParams(n_params,
                             stop_check = 1))  # Needs to be greater than 2
})

# Test for invalid stop_tol value
test_that("GetAlgoParams handles invalid stop_tol value", {
  expect_error(GetAlgoParams(n_params,
                             stop_tol = -0.1))  # Needs to be positive
})

# Test for invalid varlist type
test_that("GetAlgoParams handles invalid varlist type", {
  expect_error(GetAlgoParams(n_params,
                             varlist = 1))  # Needs to be a list
})

# Test for invalid save_rds_string format
test_that("GetAlgoParams handles invalid save_rds_string format", {
  expect_error(GetAlgoParams(n_params,
                             save_int = 1,
                             save_rds_string = 1))  # Needs to be a string
})

# Test for invalid save_int value with positive save_rds_string
test_that("GetAlgoParams handles invalid save_int with save_rds_string", {
  # Needs to be positive for saving
  expect_error(GetAlgoParams(n_params,
                             save_int = 0,
                             save_rds_string = "valid.rds"))
})
