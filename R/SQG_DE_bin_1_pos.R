#' grad_approx_fn
#'
#'
#' Computes an approximate gradient of a function using finite differences.
#'
#' @param param_indices A vector of indices specifying the parameters for which to compute the gradient.
#' @param n_diff The number of finite difference pairs to use.
#' @param current_params A matrix containing the current parameter values.
#' @param parent_indices A vector of indices specifying the parent indices for each finite difference pair.
#' @param current_weight A vector containing the current weight values.
#' @noRd
#'
#' @return A list containing:
#'   - `psi`: A normalization factor for algorithm self-scaling.
#'   - `grad_approx`: A vector of approximate gradient values.
#'
#' @details
#' This function calculates an approximate gradient by summing the differences in function values divided by the Euclidean norm of the parameter differences.
#' The normalization factor `psi` is computed to ensure proper scaling of the gradient.
#'
#' @note
#' The function includes a check to prevent division by zero and non-finite values.
#'
#' @examples
#' # NOT RUN
#' # Example usage (replace with actual data)
#' param_indices <- c(1, 2)
#' n_diff <- 10
#' current_params <- matrix(runif(20), nrow = 10, ncol = 2)
#' parent_indices <- 1:20
#' current_weight <- runif(20)
#'
#' result <- grad_approx_fn(param_indices, n_diff, current_params, parent_indices, current_weight)
#' print(result$psi)
#' print(result$grad_approx)
grad_approx_fn <-
  function(param_indices,
           n_diff,
           current_params,
           parent_indices,
           current_weight) {
    # get rid of this for loop later
    vec_diff_sum = grad_approx = numeric(length(param_indices))
    for(d in 1:n_diff){
      # difference in vector pairs
      vec_diff_temp = (current_params[parent_indices[d],param_indices] -
                         current_params[parent_indices[d+n_diff],param_indices])

      # sum up all vector differences for normalization step
      vec_diff_sum = vec_diff_sum + vec_diff_temp

      # difference in function values
      weight_diff_temp = (current_weight[parent_indices[d]] -
                            current_weight[parent_indices[d+n_diff]])

      # sum the approximate gradient vectors
      # protect against divde by zero
      # Ensure that the grad is defined
      # TODO: check the limit of this expression
      if(!(all(vec_diff_temp == 0)| weight_diff_temp == 0)){
        grad_approx = grad_approx + vec_diff_temp*(weight_diff_temp/sqrt(sum(vec_diff_temp^2)))
      }
    }

    # calculate normalization factor for algorithm self scaling
    psi_num = sqrt(sum(vec_diff_sum^2))*(1/n_diff)
    psi_den = sqrt(sum(grad_approx^2))
    psi = psi_num/psi_den
    # Ensure that the step size is defined
    # TODO: check the limit of this expression
    if(is.finite(psi_num) & psi_num == 0){
      psi = 0
    }
    # ensure that the step size down the grad in finite
    if(!is.finite(psi)){
      stop("ERROR: update is not finite!")
    }
    # return the step size and approx grad
    return(list("psi" = psi,
                "grad_approx" = grad_approx))
  }


#' SQG_DE_bin_1_pos
#'
#' @param pmem_index  Index of population (particle) member you are you are updating
#' @param current_params Current parameter values for partcle (numeric vector)
#' @param current_weight  weights for current population
#' @param objFun  function we want to minimize
#' @param n_particles number of particles
#' @param crossover_rate i.i.d. bernouli probability for updating a parameter
#' @param n_diff number of particles used to estimate the gradient
#' @param ... additional arguments for objective function
#' @noRd
#'

SQG_DE_bin_1_pos=function(pmem_index,
                          current_params,
                          params_update_ind_vec,
                          current_weight,
                          # current_like_weight = NULL,
                          # prior_function = NULL,
                          objFun,
                          scheme,
                          resample_weight,
                          step_size = .8,
                          jitter_size = 1e-6,
                          n_particles,
                          crossover_rate = 1,
                          n_diff, ... ){

  # get statistics about particle
  params_use = current_params[pmem_index,]

  # get weight to compare in greedy rule
  if(resample_weight){
    # TODO: later update to only update the prior density in blocked updating
    # if(!is.null(prior_function) & !is.null(current_like_weight)){
    #   like_weight_use = current_like_weight[pmem_index]
    #   if(all(is.finite(params_use)))weight_use = like_weight_use + prior_function(params_use, ...)
    # }else{
    # resample weight
    if(all(is.finite(params_use)))weight_use = objFun(params_use,...)
    # }
    if(!is.finite(weight_use))weight_use = Inf
    current_weight[pmem_index] <- weight_use
  }else{
    # just use the last weight
    weight_use = current_weight[pmem_index]
  }
  best_pmem_index = which.min(current_weight) # "best" specific
  # Determine how many parameters are for each particle
  len_param_use = length(params_use)

  # Determine which/how many parameters are updated in this block
  params_update = current_params[pmem_index, params_update_ind_vec]
  len_param_update = length(params_update)

  # sample parent chains
  size_to_sample <- 2*n_diff
  if(scheme == "best"){
    particles_to_sample <- c(1:n_particles)[-best_pmem_index]
  }else if (scheme == "current"){
    particles_to_sample <- c(1:n_particles)[-pmem_index]
  }else{
    particles_to_sample <- c(1:n_particles)
  }
  # sample parents
  parent_indices = sample(particles_to_sample, size=size_to_sample, replace=FALSE)

  # use binomial to sample which params to update matching crossover rate frequency
  param_idices_bool = stats::rbinom(len_param_update, prob = crossover_rate, size=1)

  # if no params selected, randomly sample 1 parameter
  if(all(param_idices_bool==0)){
    param_idices_bool[sample(x=1:len_param_update,size=1)] = 1
  }

  # indices of parameters to be updated
  param_indices = seq(1,len_param_update,by=1)[as.logical(param_idices_bool)]

  # calculate the approx grad from two parent particles
  grad_approx_list <- grad_approx_fn(param_indices,
                 n_diff,
                 current_params,
                 parent_indices,
                 current_weight)
  # get the approx grad
  grad_approx <- grad_approx_list$grad_approx
  # get the step size to go down the grad
  psi <- grad_approx_list$psi


  # Ensure that the grad and step down the gradient is finite and down the grad
  if(all(is.finite(grad_approx)) & is.finite(psi) & (psi>=0)){
    # the selected starting position of the particle before SQGDE jump
    if(scheme == "best"){
      update_index <- best_pmem_index
    }else if (scheme == "current"){
      update_index <- pmem_index
    }else {
      update_index <- parent_indices[2*n_diff+1]
    }

    # proposal: update the partcle by starting at the selected particles
    # adding a step down the approx grad
    # and adding unif noise
    params_update[param_indices] = current_params[update_index,param_indices] -
      step_size*psi*(grad_approx) + # move in the direction against the gradient
      stats::runif(length(param_indices),-jitter_size,jitter_size) # a little noise
  }
  # update the parameters in this block
  params_use[params_update_ind_vec] = params_update
  # ensure that parameters are the properly sized matrix
  params_use = matrix(params_use,1,len_param_use)

  weight_proposal = NA
  # get weight, ensure that it is finite, else assign worst value
  if(all(is.finite(params_use)))weight_proposal = objFun(params_use,...)
  if(is.na(weight_proposal))weight_proposal = Inf


  # # TODO: later update to only update the prior density in blocked updating
  # if(!is.null(prior_function) & !is.null(current_like_weight)){
  #     if(all(is.finite(params_use)))weight_proposal = objFun(params_use,...)
  #     like_weight_proposal = weight_proposal - prior_function(params_use, ...)
  # }else{
  #     if(all(is.finite(params_use)))weight_use = objFun(params_use,...)
  # }
  # if(is.na(weight_use))weight_use = Inf



  #negative greedy acceptance rule
  if(weight_proposal < weight_use) {
    current_params[pmem_index,] = params_use
    current_weight[pmem_index] = weight_proposal
    # TODO: for blocking feature
    # current_like_weight[pmem_index,] = like_weight_proposal
  }

  # TODO: for blocking feature
  # return(c(current_weight[pmem_index],
  # current_like_weight[pmem_index,],
  # current_params[pmem_index,]))
  return(c(current_weight[pmem_index],current_params[pmem_index,]))

}
