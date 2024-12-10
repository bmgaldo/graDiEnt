test_that("Purify calculates correctly", {
# Mock objective function
mock_objFun <- function(params, ...) {
  return(1)
}

# Example data
pmem_index <- 2
current_params <- matrix(runif(10), nrow = 5)
params_update_ind_vec <- c(1)
current_weight <- c(0.5, 0.8, 1.2, 0.7, 0.9)

# Run the test
# debugonce(Purify)
updated_values <- Purify(pmem_index,
                         current_params,
                         params_update_ind_vec,
                         current_weight, mock_objFun)

# Check if weight and params are updated
expect_equal(updated_values[1], 1)  # New weight
# parmetes do not change
expect_true(all.equal(updated_values[-1], current_params[2, ]))
})
