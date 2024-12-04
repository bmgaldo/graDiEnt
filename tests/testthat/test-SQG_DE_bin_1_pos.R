library(testthat)

test_that("resample_weight updates weight if resample_weight is TRUE", {
  # Mock objective function
  mock_obj_fun <- function(params, ...) {
    return(sum(params^2))
  }

  # Set up test data that produces values where a proposal is no different than
  # the former candidate
  current_params <- matrix(rep(c(1, 2), each = 4), nrow = 4)
  current_weight <- rep(10, 4)
  pmem_index <- 1
  resample_weight <- TRUE

  # Call the function
  updated_weight_params <- SQG_DE_bin_1_pos(pmem_index = pmem_index,
                                     current_params = current_params,
                                     params_update_ind_vec = c(1, 2),
                                     current_weight = current_weight,
                                     objFun = mock_obj_fun,
                                     scheme = "best",
                                     n_diff = 1,
                                     n_particles = 4,
                                     resample_weight = resample_weight,
                                     jitter_size = 0)

  # Expected weight after objective function call
  expected_weight <- mock_obj_fun(current_params[pmem_index,])

  # Check if weight is updated and matches expected value
  expect_equal(updated_weight_params[1], expected_weight)
})

test_that("resample_weight uses current weight if resample_weight is FALSE", {
  # Mock objective function (not called in this test)
  mock_obj_fun <- function(params, ...) {
    return(sum(params^2))
  }

  # Set up test data that produces values where a proposal is no different than
  # the former candidate
  current_params <- matrix(rep(c(1, 2), each = 4), nrow = 4)
  current_weight <- rep(1, 4)
  pmem_index <- 1
  resample_weight <- FALSE

  # Call the function
  updated_weight_params <- SQG_DE_bin_1_pos(pmem_index = pmem_index,
                                     current_params = current_params,
                                     params_update_ind_vec = c(1, 2),
                                     current_weight = current_weight,
                                     objFun = mock_obj_fun,
                                     scheme = "best",
                                     n_diff = 1,
                                     n_particles = 4,
                                     resample_weight = resample_weight)

  # Check if weight is not updated and remains the same
  expect_equal(updated_weight_params[1], current_weight[1])
})

test_that("parent chain sampling is correct", {

  # Mock objective function (not called in this test)
  mock_obj_fun <- function(params, ...) {
    return(sum(params^2))
  }

  # Set up test data that produces values where a proposal is no different than
  # the former candidate and there is a best particle
  current_params <- matrix(rep(c(1, 2), each = 4), nrow = 4)
  current_weight <- rep(1000, 4)
  current_weight[4] <- 0 # the best index is 4
  current_params[4, ] <- c(0, 0)
  pmem_index <- 1
  resample_weight <- FALSE

  # Test "best" scheme
  scheme <- "best"
  # Test should update the index 1 to the parameterization of 4
  updated_weight_params <- SQG_DE_bin_1_pos(pmem_index = pmem_index,
                                     current_params = current_params,
                                     params_update_ind_vec = c(1, 2),
                                     current_weight = current_weight,
                                     objFun = mock_obj_fun,
                                     scheme = scheme,
                                     n_diff = 1,
                                     n_particles = 4,
                                     resample_weight = resample_weight,
                                     jitter_size = 0)
  expect_equal(updated_weight_params, c(current_weight[4], current_params[4, ]))


  # Set up test data that produces values where a proposal is no different than
  # the former candidate
  current_params <- matrix(rep(c(1, 2), each = 4), nrow = 4)
  current_weight <- rep(1000, 4)
  current_params[1, ] <- c(10, 10)

  # Test "current" scheme
  scheme <- "current"
  updated_weight_params <- SQG_DE_bin_1_pos(pmem_index = pmem_index,
                                            current_params = current_params,
                                            params_update_ind_vec = c(1, 2),
                                            current_weight = current_weight,
                                            objFun = mock_obj_fun,
                                            scheme = scheme,
                                            n_diff = 1,
                                            n_particles = 4,
                                            resample_weight = resample_weight,
                                            jitter_size = 0)
  expect_equal(updated_weight_params, c(mock_obj_fun(current_params[1, ]), current_params[1, ]))

  # Test "rand" scheme
  scheme <- "rand"
  # IDK how this can be tested for now
  # updated_weight_params <- SQG_DE_bin_1_pos(pmem_index = pmem_index,
  #                                           current_params = current_params,
  #                                           params_update_ind_vec = c(1, 2),
  #                                           current_weight = current_weight,
  #                                           objFun = mock_obj_fun,
  #                                           scheme = scheme,
  #                                           n_diff = 1,
  #                                           n_particles = 4,
  #                                           resample_weight = resample_weight,
  #                                           jitter_size = 0)
  # expect_equal(updated_weight_params, c(mock_obj_fun(current_params[1, ]), current_params[1, ]))


})

test_that("parameter selection is correct", {

  # Mock objective function (not called in this test)
  mock_obj_fun <- function(params, ...) {
    return(sum(params^2))
  }

  # Set up test data that produces values where a proposal is no different than
  # the former candidate
  num_params <- 1e2
  current_params <- matrix(rnorm(num_params*4), nrow = 4)
  current_weight <- 100000 + 0:3
  pmem_index <- 1
  resample_weight <- FALSE
  crossover_rate <- 0

  scheme <- "current"
  # debugonce(SQG_DE_bin_1_pos)
  updated_weight_params <- SQG_DE_bin_1_pos(pmem_index = pmem_index,
                                            current_params = current_params,
                                            params_update_ind_vec = 1:num_params,
                                            current_weight = current_weight,
                                            objFun = mock_obj_fun,
                                            scheme = scheme,
                                            n_diff = 1,
                                            n_particles = 4,
                                            resample_weight = resample_weight,
                                            jitter_size = 0,
                                            crossover_rate = crossover_rate)
  expect_lt(updated_weight_params[1], current_weight[1])
  expect_equal(sum(updated_weight_params[-1] == current_params[1,]),  num_params-1)


  num_params <- 1e7
  current_params <- matrix(rnorm(num_params*4), nrow = 4)
  current_weight <- 1e9 + 0:3
  pmem_index <- 1
  resample_weight <- FALSE
  crossover_rate <- .5
  updated_weight_params <- SQG_DE_bin_1_pos(pmem_index = pmem_index,
                                            current_params = current_params,
                                            params_update_ind_vec = 1:num_params,
                                            current_weight = current_weight,
                                            objFun = mock_obj_fun,
                                            scheme = scheme,
                                            n_diff = 1,
                                            n_particles = 4,
                                            resample_weight = resample_weight,
                                            jitter_size = 0,
                                            crossover_rate = crossover_rate)
  expect_lt(updated_weight_params[1], current_weight[1])
  diff <-
    sum(updated_weight_params[-1] == current_params[1,]) / num_params -
    crossover_rate
  expect_lt(abs(diff), .01)
})

library(testthat)

test_that("grad_approx_fn calculates correctly", {
  # Set up test data
  current_params <- matrix(1:16, nrow = 4, ncol = 4)
  current_weight <- c(1, 2, 3, 4)
  parent_indices <- c(1, 2, 3, 4)
  n_diff <- 2
  param_indices <- c(1, 3)

  # Calculate gradient approximation using the function
  result <- grad_approx_fn(param_indices, n_diff, current_params, parent_indices, current_weight)

  # Calculate gradient approximation manually
  vec_diff_sum <- numeric(length(param_indices))
  grad_approx_manual <- numeric(length(param_indices))
  for(d in 1:n_diff){
    vec_diff_temp <- current_params[parent_indices[d], param_indices] - current_params[parent_indices[d + n_diff], param_indices]
    vec_diff_sum <- vec_diff_sum + vec_diff_temp
    weight_diff_temp <- current_weight[parent_indices[d]] - current_weight[parent_indices[d + n_diff]]
    if(!(all(vec_diff_temp == 0) | weight_diff_temp == 0)) {
      grad_approx_manual <- grad_approx_manual + vec_diff_temp * (weight_diff_temp / sqrt(sum(vec_diff_temp^2)))
    }
  }

  # Calculate psi manually
  psi_num <- sqrt(sum(vec_diff_sum^2)) * (1 / n_diff)
  psi_den <- sqrt(sum(grad_approx_manual^2))
  psi_manual <- psi_num / psi_den

  # Compare the two results
  expect_equal(result$grad_approx, grad_approx_manual, tolerance = 1e-10)
  expect_equal(result$psi, psi_manual, tolerance = 1e-10)
})

