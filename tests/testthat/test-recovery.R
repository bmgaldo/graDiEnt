library(graDiEnt)

obj <- function(x) sum((x - 1)^2)

# ── GetAlgoParams validation ──────────────────────────────────────────────────

test_that("GetAlgoParams accepts recovery_path = NULL", {
  cp <- GetAlgoParams(n_params = 2, recovery_path = NULL)
  expect_null(cp$recovery_path)
})

test_that("GetAlgoParams accepts a character scalar recovery_path", {
  cp <- GetAlgoParams(n_params = 2, recovery_path = tempfile(fileext = ".rds"))
  expect_true(is.character(cp$recovery_path))
})

test_that("GetAlgoParams rejects non-character recovery_path", {
  expect_error(
    GetAlgoParams(n_params = 2, recovery_path = 123),
    "character scalar"
  )
})

test_that("GetAlgoParams rejects non-scalar recovery_path", {
  expect_error(
    GetAlgoParams(n_params = 2, recovery_path = c("a.rds", "b.rds")),
    "character scalar"
  )
})

test_that("GetAlgoParams rejects recovery_path without .rds extension", {
  expect_error(
    GetAlgoParams(n_params = 2, recovery_path = "output.txt"),
    "\\.rds"
  )
})

test_that("GetAlgoParams rejects recovery_path with no extension", {
  expect_error(
    GetAlgoParams(n_params = 2, recovery_path = "myfile"),
    "\\.rds"
  )
})

test_that("GetAlgoParams accepts .RDS extension (case-insensitive)", {
  expect_no_error(GetAlgoParams(n_params = 2, recovery_path = "output.RDS"))
})

# ── recovery_freq validation ──────────────────────────────────────────────────

test_that("GetAlgoParams accepts recovery_freq = 1", {
  cp <- GetAlgoParams(n_params = 2, recovery_freq = 1)
  expect_equal(cp$recovery_freq, 1L)
})

test_that("GetAlgoParams accepts recovery_freq = 10", {
  cp <- GetAlgoParams(n_params = 2, recovery_freq = 10)
  expect_equal(cp$recovery_freq, 10L)
})

test_that("GetAlgoParams rejects recovery_freq = 0", {
  expect_error(GetAlgoParams(n_params = 2, recovery_freq = 0), "positive integer")
})

test_that("GetAlgoParams rejects negative recovery_freq", {
  expect_error(GetAlgoParams(n_params = 2, recovery_freq = -5), "positive integer")
})

test_that("GetAlgoParams rejects non-scalar recovery_freq", {
  expect_error(GetAlgoParams(n_params = 2, recovery_freq = c(1, 2)), "positive integer")
})

# ── recovery file written during run ─────────────────────────────────────────

test_that("recovery file is created and readable after run", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 1, n_iter = 50, n_particles = 6,
                              n_diff = 1, recovery_path = f))
  )

  expect_true(file.exists(f))
  rec <- readRDS(f)
  expect_true(is.list(rec))
})

test_that("recovery file has correct structure", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 8,
                              n_diff = 1, recovery_path = f))
  )

  rec <- readRDS(f)
  expect_named(rec, c("solution", "weight", "last_particles", "last_weights",
                      "converged", "iter_completed"),
               ignore.order = TRUE)
})

test_that("recovery file solution and weight are finite", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 2, n_iter = 50, n_particles = 8,
                              n_diff = 1, recovery_path = f))
  )

  rec <- readRDS(f)
  expect_true(all(is.finite(rec$solution)))
  expect_true(is.finite(rec$weight))
})

test_that("recovery last_particles has correct dimensions", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 3, n_iter = 50, n_particles = 12,
                              n_diff = 1, recovery_path = f))
  )

  rec <- readRDS(f)
  expect_equal(dim(rec$last_particles), c(12L, 3L))
  expect_equal(length(rec$last_weights), 12L)
})

test_that("recovery iter_completed equals n_iter when run completes", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 1, n_iter = 30, n_particles = 6,
                              n_diff = 1, recovery_path = f,
                              stop_check = 40))  # prevent early stop
  )

  rec <- readRDS(f)
  expect_equal(rec$iter_completed, 30L)
})

test_that("recovery converged field is FALSE", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 1, n_iter = 30, n_particles = 6,
                              n_diff = 1, recovery_path = f,
                              stop_check = 40))
  )

  rec <- readRDS(f)
  expect_false(rec$converged)
})

# ── recovery file does not contain trace arrays ───────────────────────────────

test_that("recovery file excludes trace arrays even when return_trace = TRUE", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 1, n_iter = 30, n_particles = 6,
                              n_diff = 1, recovery_path = f,
                              return_trace = TRUE, stop_check = 40))
  )

  rec <- readRDS(f)
  expect_false("particles_trace" %in% names(rec))
  expect_false("weights_trace"   %in% names(rec))
})

# ── recovery_freq controls save timing ───────────────────────────────────────

test_that("recovery file not written when n_iter < recovery_freq", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 1, n_iter = 5, n_particles = 6,
                              n_diff = 1, recovery_path = f, recovery_freq = 10,
                              stop_check = 20))
  )
  expect_false(file.exists(f))
})

test_that("recovery file written when iter reaches recovery_freq multiple", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))

  set.seed(1)
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 1, n_iter = 20, n_particles = 6,
                              n_diff = 1, recovery_path = f, recovery_freq = 10,
                              stop_check = 30))
  )
  expect_true(file.exists(f))
  rec <- readRDS(f)
  expect_true(rec$iter_completed %% 10 == 0)
})

# ── no file written when recovery_path = NULL ─────────────────────────────────

test_that("no rds file written when recovery_path is NULL", {
  set.seed(1)
  before <- list.files(tempdir(), pattern = "\\.rds$")
  suppressMessages(
    optim_SQGDE(obj,
                GetAlgoParams(n_params = 1, n_iter = 20, n_particles = 6,
                              n_diff = 1, recovery_path = NULL))
  )
  after <- list.files(tempdir(), pattern = "\\.rds$")
  expect_equal(before, after)
})
