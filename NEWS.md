# graDiEnt 1.1.0

## New features

* **Box constraints**: `GetAlgoParams` gains `lower` and `upper` parameters for per-parameter bounds. Enforced via `bounds_type = "reflect"` (default, mirrors proposals back across boundary) or `"clip"` (truncates to boundary).
* **Warm start**: `optim_SQGDE` gains `warm_start` parameter. Pass the output of a previous run to seed the population and continue optimization without re-initializing.
* **Disk recovery**: `GetAlgoParams` gains `recovery_path` (`.rds` file path) and `recovery_freq` (save interval in iterations). Partial results are written to disk periodically via `saveRDS`, surviving interrupts and OOM crashes. Load with `readRDS(recovery_path)`.
* **Progress message control**: `GetAlgoParams` gains `trace_print_freq`. Set to an integer to print progress every N iterations, or `Inf` to suppress all iteration messages.
* **Unified adaptation schemes**: `rand`, `current`, and `best` DE schemes are now selected via `adapt_scheme` in `GetAlgoParams` rather than separate internal functions.
* **`crossover_rate = 0`**: Now accepted; selects exactly one random parameter to update per iteration.
* **`jitter_size = 0`**: Now accepted; disables jitter noise.
* **Block coordinate updating**: `GetAlgoParams` gains `param_block_list`, an optional list of integer vectors partitioning `1:n_params` into blocks. The optimizer cycles through blocks in order across iterations, updating only the current block's parameters instead of using random crossover selection. Useful for high-dimensional problems where updating parameter subsets per iteration is more efficient.

## Bug fixes

* Fixed division-by-zero in `grad_approx_fn` when two particles are identical or have equal weights (zero-diff guard).
* Fixed `rand` scheme incorrectly sampling `2*n_diff` parents instead of `2*n_diff + 1`, causing out-of-bounds index access.
* Fixed `runif` jitter applied over `n_params` elements instead of `length(param_indices)` when `crossover_rate < 1`.
* Fixed `GetAlgoParams` validation for `n_params`, `n_particles`, and `n_iter`: finiteness is now checked before `as.integer()` coercion, preventing silent `NA` production and spurious "condition has length > 1" errors on vector inputs.

# graDiEnt 1.0.1
* Package uses "message()" instead of "print()" to communicate with user.
* minor fix to documentation.

# graDiEnt 1.0.0
* First public release of package.
