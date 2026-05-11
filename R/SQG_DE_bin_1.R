# Internal helper: stochastic quasi-gradient approximation.
# Returns psi (self-scaling factor) and grad_approx (approximate gradient vector).
# Skips any pair where the parameter difference or weight difference is zero
# to avoid division by zero.
grad_approx_fn = function(param_indices, n_diff, current_params,
                           parent_indices, current_weight) {
  vec_diff_sum = grad_approx = numeric(length(param_indices))
  for (d in 1:n_diff) {
    vec_diff_temp = (current_params[parent_indices[d],     param_indices] -
                     current_params[parent_indices[d+n_diff], param_indices])
    vec_diff_sum  = vec_diff_sum + vec_diff_temp
    weight_diff_temp = (current_weight[parent_indices[d]] -
                        current_weight[parent_indices[d+n_diff]])
    if (!(all(vec_diff_temp == 0) | weight_diff_temp == 0)) {
      grad_approx = grad_approx +
        vec_diff_temp * (weight_diff_temp / sqrt(sum(vec_diff_temp^2)))
    }
  }
  psi_num = sqrt(sum(vec_diff_sum^2)) * (1/n_diff)
  psi_den = sqrt(sum(grad_approx^2))
  list(psi = psi_num / psi_den, grad_approx = grad_approx)
}


#' SQG_DE_bin_1
#'
#' @param pmem_index Index of population member being updated.
#' @param current_params Current parameter matrix (n_particles x n_params).
#' @param current_weight Current weight vector (length n_particles).
#' @param objFun Objective function to minimize.
#' @param scheme Adaptation scheme: 'rand', 'best', or 'current'.
#' @param step_size Scaling factor for gradient step.
#' @param jitter_size Half-width of uniform jitter added to proposals.
#' @param n_particles Number of particles.
#' @param crossover_rate Bernoulli probability of updating each parameter.
#' @param lower Lower bound vector (length n_params).
#' @param upper Upper bound vector (length n_params).
#' @param bounds_type Bounds enforcement: 'reflect' or 'clip'.
#' @param n_diff Number of particle pairs used for gradient approximation.
#' @param ... Additional arguments passed to objFun.
#' @noRd

SQG_DE_bin_1 = function(pmem_index,
                         current_params,
                         current_weight,
                         objFun,
                         scheme       = 'rand',
                         step_size    = .8,
                         jitter_size  = 1e-6,
                         n_particles,
                         crossover_rate = 1,
                         lower        = -Inf,
                         upper        =  Inf,
                         bounds_type  = 'reflect',
                         n_diff, ...) {

  weight_use       = current_weight[pmem_index]
  params_use       = current_params[pmem_index, ]
  len_param_use    = length(params_use)
  best_pmem_index  = which.min(current_weight)

  # scheme determines: parent pool, sample size, and base position for proposal
  if (scheme == 'best') {
    parent_indices = sample(c(1:n_particles)[-best_pmem_index],
                            size = 2*n_diff, replace = FALSE)
    base_index = best_pmem_index
  } else if (scheme == 'current') {
    parent_indices = sample(c(1:n_particles)[-pmem_index],
                            size = 2*n_diff, replace = FALSE)
    base_index = pmem_index
  } else {
    parent_indices = sample(1:n_particles, size = 2*n_diff+1, replace = FALSE)
    base_index = parent_indices[2*n_diff+1]
  }

  # crossover: select which parameters to update
  param_idices_bool = stats::rbinom(len_param_use, prob = crossover_rate, size = 1)
  if (all(param_idices_bool == 0)) {
    param_idices_bool[sample(x = 1:len_param_use, size = 1)] = 1
  }
  param_indices = seq(1, len_param_use, by = 1)[as.logical(param_idices_bool)]

  # approximate gradient and self-scaling factor
  ga           = grad_approx_fn(param_indices, n_diff, current_params,
                                parent_indices, current_weight)
  grad_approx  = ga$grad_approx
  psi          = ga$psi

  if (all(is.finite(grad_approx)) & is.finite(psi) & (psi > 0)) {
    params_use[param_indices] = current_params[base_index, param_indices] -
      step_size * psi * grad_approx +
      stats::runif(length(param_indices), -jitter_size, jitter_size)
  }

  params_use = apply_bounds(params_use, lower, upper, bounds_type)
  params_use = matrix(params_use, 1, len_param_use)

  weight_proposal = NA
  if (all(is.finite(params_use))) weight_proposal = objFun(params_use, ...)
  if (is.na(weight_proposal)) weight_proposal = Inf

  if (weight_proposal < weight_use) {
    current_params[pmem_index, ] = params_use
    current_weight[pmem_index]   = weight_proposal
  }

  return(c(current_weight[pmem_index], current_params[pmem_index, ]))
}
