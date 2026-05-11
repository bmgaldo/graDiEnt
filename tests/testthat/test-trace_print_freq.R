library(graDiEnt)

obj <- function(x) sum((x - 1)^2)

# helper: run optimizer and return only "iter N/M" progress messages
capture_iter_messages <- function(control_params) {
  msgs <- character(0)
  withCallingHandlers(
    optim_SQGDE(obj, control_params),
    message = function(m) {
      txt <- conditionMessage(m)
      if (grepl("^iter ", txt)) msgs <<- c(msgs, txt)
      invokeRestart("muffleMessage")
    }
  )
  msgs
}

# ── Inf suppresses all iter messages ─────────────────────────────────────────

test_that("trace_print_freq = Inf produces no iter progress messages", {
  set.seed(1)
  msgs <- capture_iter_messages(
    GetAlgoParams(n_params = 1, n_iter = 200, n_particles = 6,
                  n_diff = 1, trace_print_freq = Inf)
  )
  expect_length(msgs, 0)
})

# ── custom freq prints at correct multiples ───────────────────────────────────

test_that("trace_print_freq = 10 prints at iterations 10 and 20", {
  set.seed(1)
  msgs <- capture_iter_messages(
    GetAlgoParams(n_params = 1, n_iter = 25, n_particles = 6,
                  n_diff = 1, trace_print_freq = 10, stop_check = 30)
  )
  expect_length(msgs, 2)
  expect_match(msgs[1], "^iter 10/")
  expect_match(msgs[2], "^iter 20/")
})

test_that("trace_print_freq = 50 prints at iterations 50 and 100", {
  set.seed(1)
  msgs <- capture_iter_messages(
    GetAlgoParams(n_params = 1, n_iter = 110, n_particles = 6,
                  n_diff = 1, trace_print_freq = 50, stop_check = 120)
  )
  expect_length(msgs, 2)
  expect_match(msgs[1], "^iter 50/")
  expect_match(msgs[2], "^iter 100/")
})

# ── validation errors ─────────────────────────────────────────────────────────

test_that("GetAlgoParams rejects trace_print_freq = 0", {
  expect_error(
    GetAlgoParams(n_params = 2, trace_print_freq = 0),
    "positive integer or Inf"
  )
})

test_that("GetAlgoParams rejects negative trace_print_freq", {
  expect_error(
    GetAlgoParams(n_params = 2, trace_print_freq = -10),
    "positive integer or Inf"
  )
})

test_that("GetAlgoParams rejects non-scalar trace_print_freq", {
  expect_error(
    GetAlgoParams(n_params = 2, trace_print_freq = c(10, 20)),
    "scalar"
  )
})
