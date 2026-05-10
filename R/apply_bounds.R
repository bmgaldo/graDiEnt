# Internal helper: enforce box constraints on a parameter vector.
# "clip"    — truncate to [lower, upper]
# "reflect" — reflect off each boundary once, then clip (handles large steps)
apply_bounds <- function(x, lower, upper, type) {
  if (type == "reflect") {
    below <- is.finite(lower) & x < lower
    x[below] <- 2 * lower[below] - x[below]
    above <- is.finite(upper) & x > upper
    x[above] <- 2 * upper[above] - x[above]
  }
  pmax(lower, pmin(upper, x))
}
