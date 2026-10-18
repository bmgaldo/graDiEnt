## Test environments

* local: macOS arm64 (Apple Silicon), R-release 4.4.1 — 0 errors, 0 warnings
* GitHub Actions: ubuntu-latest, R-release — OK
* R-hub v2 (GitHub Actions):
  * ubuntu-clang (R-devel) — Status: OK
  * ubuntu-gcc12 (R-devel) — Status: OK
  * ubuntu-next (R 4.6.1 Patched) — Status: OK
  * ubuntu-release (R 4.6.1) — Status: OK
  * macos-arm64 (R-devel) — Status: OK
* win-builder: r-devel, r-release, r-oldrel — 0 errors, 0 warnings, 0 notes

## R CMD check results

0 errors | 0 warnings | 0 notes

If the incoming check flags possibly misspelled words in DESCRIPTION, these are
the expected candidates:

* Baldanzini, Sala, Pierini — last names of authors of the cited paper
  (Sala, Baldanzini, and Pierini 2018, doi:10.1007/978-3-319-72926-8_27)
* SQG — abbreviation for the algorithm name (Stochastic Quasi-Gradient)
* optim — name of a function in base R whose interface convention this package follows

These words are spelled correctly or are proper nouns.

## Changes since last CRAN release (1.0.1, 2022-05-10)

This is a major feature update to 1.1.0:

* New: box constraints (`lower`/`upper` with `reflect` or `clip` enforcement),
  warm start, disk recovery (`recovery_path`/`recovery_freq`), progress message
  control (`trace_print_freq`), unified DE adaptation schemes (`adapt_scheme`),
  block coordinate updating (`param_block_list`), `crossover_rate = 0` and
  `jitter_size = 0` support.
* Fixed: bugs in gradient approximation (division-by-zero guard), parent
  sampling under the `rand` scheme, jitter dimensions with `crossover_rate < 1`,
  parameter validation, and an early-stopping crash when `thin > 1`
  (see NEWS.md for the full list).
* Test suite expanded from 4 tests to 228 tests covering all exported and
  internal functions.