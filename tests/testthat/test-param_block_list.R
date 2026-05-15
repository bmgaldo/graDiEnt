library(graDiEnt)

obj <- function(x) sum((x - 1)^2)

# ── GetAlgoParams validation ──────────────────────────────────────────────────

test_that("param_block_list = NULL accepted (default)", {
  cp <- GetAlgoParams(n_params = 4)
  expect_null(cp$param_block_list)
})

test_that("valid partition accepted", {
  expect_no_error(
    GetAlgoParams(n_params = 4, param_block_list = list(1:2, 3:4))
  )
})

test_that("single-block partition accepted", {
  expect_no_error(
    GetAlgoParams(n_params = 3, param_block_list = list(1:3))
  )
})

test_that("unequal block sizes accepted", {
  expect_no_error(
    GetAlgoParams(n_params = 4, param_block_list = list(c(1L, 2L, 3L), 4L))
  )
})

test_that("indices coerced to integer", {
  cp <- GetAlgoParams(n_params = 4, param_block_list = list(c(1.0, 2.0), c(3.0, 4.0)))
  expect_true(is.integer(cp$param_block_list[[1]]))
})

test_that("non-list param_block_list errors", {
  expect_error(
    GetAlgoParams(n_params = 4, param_block_list = 1:4),
    "list"
  )
})

test_that("index out of range errors", {
  expect_error(
    GetAlgoParams(n_params = 4, param_block_list = list(1:3, c(4L, 5L))),
    "1:n_params"
  )
})

test_that("duplicate index errors", {
  expect_error(
    GetAlgoParams(n_params = 4, param_block_list = list(c(1L, 2L), c(2L, 3L, 4L))),
    "exactly once"
  )
})

test_that("missing index errors", {
  expect_error(
    GetAlgoParams(n_params = 4, param_block_list = list(1L, 2L, 3L)),
    "exactly once"
  )
})

test_that("empty list errors", {
  expect_error(
    GetAlgoParams(n_params = 4, param_block_list = list()),
    "non-empty"
  )
})

# ── optimizer runs with param_block_list ──────────────────────────────────────

test_that("optimizer runs without error with 2-block partition", {
  set.seed(1)
  expect_no_error(
    suppressMessages(
      optim_SQGDE(obj,
                  GetAlgoParams(n_params = 4, n_iter = 100, n_particles = 12,
                                n_diff = 2,
                                param_block_list = list(1:2, 3:4)))
    )
  )
})

test_that("optimizer runs without error with single-element blocks", {
  set.seed(1)
  expect_no_error(
    suppressMessages(
      optim_SQGDE(obj,
                  GetAlgoParams(n_params = 3, n_iter = 100, n_particles = 12,
                                n_diff = 2,
                                param_block_list = list(1L, 2L, 3L)))
    )
  )
})

test_that("optimizer with param_block_list returns finite solution", {
  set.seed(1)
  out <- suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 4, n_iter = 200, n_particles = 12,
                              n_diff = 2, init_sd = 1,
                              param_block_list = list(1:2, 3:4)))
  )
  expect_true(all(is.finite(out$solution)))
  expect_true(is.finite(out$weight))
})

test_that("block cycling: n_params=1 with single block works", {
  set.seed(1)
  out <- suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 1, n_iter = 50, n_particles = 6,
                              n_diff = 1,
                              param_block_list = list(1L)))
  )
  expect_true(is.finite(out$solution))
})
