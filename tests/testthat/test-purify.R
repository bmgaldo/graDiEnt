library(testthat)

# Define a simple objective function for testing
test_obj_fun <- function(x) {
  return(sum((x - 3)^2))
}

test_that("Purify updates weights and parameters correctly", {

  # Initialize test data
  current_params <- matrix(c(2.5, 3.5,
                             4.5, 5.5), nrow = 2, byrow = TRUE)
  current_weight <- c(10, 20)
  n_particles <- 2

  # Test Purify function
  result <- Purify(
    pmem_index = 1,
    current_params = current_params,
    current_weight = current_weight,
    objFun = test_obj_fun,
    n_particles = n_particles
  )

  # Expected weight after running objFun on current_params[1,]
  expected_weight <- test_obj_fun(current_params[1, ])

  # Check that the weight and parameters were updated correctly
  expect_equal(result[1], expected_weight)
  expect_equal(result[-1], current_params[1, ])
})

test_that("Purify handles NA and Inf values correctly", {

  # Test case where objFun returns NA
  na_obj_fun <- function(x) NA
  current_params <- matrix(c(2.5, 3.5), nrow = 1)
  current_weight <- c(10)
  n_particles <- 1

  result <- Purify(
    pmem_index = 1,
    current_params = current_params,
    current_weight = current_weight,
    objFun = na_obj_fun,
    n_particles = n_particles
  )

  # Expect weight to remain Inf when objFun returns NA
  expect_equal(result[1], 10)

  # Test case where objFun returns Inf
  inf_obj_fun <- function(x) Inf
  result <- Purify(
    pmem_index = 1,
    current_params = current_params,
    current_weight = current_weight,
    objFun = inf_obj_fun,
    n_particles = n_particles
  )

  # Expect weight to remain Inf when objFun returns Inf
  expect_equal(result[1], 10)
})

test_that("Purify returns correct output structure", {

  # Initialize test data
  current_params <- matrix(c(2.5, 3.5), nrow = 1)
  current_weight <- c(10)
  n_particles <- 1

  # Test Purify function
  result <- Purify(
    pmem_index = 1,
    current_params = current_params,
    current_weight = current_weight,
    objFun = test_obj_fun,
    n_particles = n_particles
  )

  # Check that the result has the correct length
  expect_equal(length(result), ncol(current_params) + 1)

  # Check that the first element is a numeric weight
  expect_type(result[1], "double")

  # Check that the remaining elements are numeric parameters
  expect_type(result[-1], "double")
})
