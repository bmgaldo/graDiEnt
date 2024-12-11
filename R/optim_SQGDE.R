#' optim_SQGDE
#'
#' @description Uses Stochastic Quasi-Gradient Differential Evolution (SQG-DE;
#' Sala, Baldanzini, and Pierini, 2018) to minimize a specified objective function `f(x)`.
#' To maximize a function, pass `g(x) = -f(x)` as the `ObjFun` argument.
#'
#' @param ObjFun A scalar-returning function or a list of functions to minimize. Each function
#' should accept a real-valued vector of length `n_params` as input. When providing a list of functions,
#' block updating is enabled, iterating through each function to update only a subset of parameters conditional.
#'
#' @param control_params A list of control parameters for the SQG-DE algorithm. The primary required
#' parameter is `n_params`. For more details, see the documentation of \code{\link{GetAlgoParams}}.
#'
#' @param ... Additional arguments to be passed to `ObjFun`.
#'
#' @return A list containing:
#' - `solution`: The best solution found by the algorithm.
#' - `weight`: The corresponding weight (i.e., `f(solution)`).
#' - `particles_trace` (optional): Trace of particle positions over iterations (if `return_trace = TRUE`).
#' - `weights_trace` (optional): Trace of weights over iterations (if `return_trace = TRUE`).
#' - `converged`: Logical, indicating whether the convergence criterion was met.
#'
#' @export
#' @md
#' @examples
#' ##############
#' # Maximum Likelihood Example
#' ##############
#'
#' # Simulate data
#' dataExample <- matrix(rnorm(1000, c(-1, 1, 0, 1), c(1, 1, 1, 1)), ncol = 4, byrow = TRUE)
#'
#' # Define parameter names
#' param_names_example <- c("mu_1", "mu_2", "mu_3", "mu_4")
#'
#' # Define negative log-likelihood function
#' ExampleObjFun <- function(x, data, param_names) {
#'   names(x) <- param_names
#'   -sum(dnorm(data[, 1], x["mu_1"], sd = 1, log = TRUE) +
#'        dnorm(data[, 2], x["mu_2"], sd = 1, log = TRUE) +
#'        dnorm(data[, 3], x["mu_3"], sd = 1, log = TRUE) +
#'        dnorm(data[, 4], x["mu_4"], sd = 1, log = TRUE))
#' }
#'
#' # Run optimization
#' out <- optim_SQGDE(
#'   ObjFun = ExampleObjFun,
#'   control_params = GetAlgoParams(
#'     n_params = length(param_names_example),
#'     n_iter = 250,
#'     n_particles = 12,
#'     n_diff = 2,
#'     return_trace = TRUE
#'   ),
#'   data = dataExample,
#'   param_names = param_names_example
#' )
#'
#' # Plot particle trajectories
#' par(mfrow = c(2, 2))
#' matplot(out$particles_trace[, , 1], type = 'l')
#' matplot(out$particles_trace[, , 2], type = 'l')
#' matplot(out$particles_trace[, , 3], type = 'l')
#' matplot(out$particles_trace[, , 4], type = 'l')
#'
#' # SQG-DE solution
#' out$solution
#'
#' # Analytical solution
#' apply(dataExample, 2, mean)


optim_SQGDE = function(ObjFun,
                       # prior_function_list = NULL,
                       control_params = GetAlgoParams(), ...){


  # ObjFun
  ### catch errors
  if(length(ObjFun) == 1){
    if(!is.function(ObjFun)){
      stop('ERROR: ObjFun is not A FUNCTION or list of functions!')
    }
    ObjFun <- list(ObjFun)
  }else{
    if(!is.list(ObjFun)){
      stop('ERROR: ObjFun is not a function or LIST of functions!')
    }
    for(f in seq_along(ObjFun)){
      if(!is.function(ObjFun[[f]])){
        stop(paste0('ERROR: The ', f, 'element of ObjFun is not a function!'))
      }
    }
  }
  if(length(ObjFun) != length(control_params$param_ind_to_update_list)){
    stop(paste0('ERROR: The element of ObjFun has a length of ',
                length(ObjFun),
                ' which is not the same length as the param_ind_to_update_list, ',
                length(control_params$param_ind_to_update_list), '!'))
  }

  # prior_function_list
  ### catch errors
  # if(!is.null(prior_function_list)){
  #   if(length(prior_function_list) == 1){
  #     if(!is.function(prior_function_list)){
  #       stop('ERROR: prior_function_list is not A FUNCTION or list of functions!')
  #     }
  #     prior_function_list <- list(prior_function_list)
  #   }else{
  #     if(!is.list(prior_function_list)){
  #       stop('ERROR: prior_function_list is not a function or LIST of functions!')
  #     }
  #     for(f in seq_along(prior_function_list)){
  #       if(!is.function(prior_function_list[[f]])){
  #         stop(paste0('ERROR: The ', f, 'element of prior_function_list is not a function!'))
  #       }
  #     }
  #   }
  #   if(length(prior_function_list) != length(ObjFun)){
  #     stop(paste0('ERROR: The prior_function_list is not the same length as ObjFun!'))
  #   }
  # }

  # create memory structures for storing particle trajectories
  particles = array(NA,
                    dim = c(control_params$n_iters_per_particle,
                            control_params$n_particles,
                            control_params$n_params))
  weights = array(NA,
                  dim = c(control_params$n_iters_per_particle,
                          control_params$n_particles,
                          ObjFun |> length()))
  # TODO: for blocking
  # like_weights = array(NA,
  #                      dim = c(control_params$n_iters_per_particle,
  #                              control_params$n_particles,
  #                              ObjFun |> length()))

  # cluster initialization
  if(!control_params$parallel_type=='none'){
    message(paste0("initalizing ",
                   control_params$parallel_type, " cluser with ",
                   control_params$n_cores_use, " cores"))

    doParallel::registerDoParallel(control_params$n_cores_use)
    if(control_params$parallel_type == "PSOCK"){
      cl_use =
        parallel::makePSOCKcluster(
          names = control_params$n_cores_use,
          outfile = control_params$outfile)
    }else if(control_params$parallel_type == "FORK"){
      cl_use =
        parallel::makeForkCluster(control_params$n_cores_use,
                                  outfile = control_params$outfile)
    }

    if(!is.null(control_params$parallel_seed)){
      parallel::clusterSetRNGStream(cl_use, control_params$parallel_seed)
    }
    parallel::clusterExport(cl_use,
                            varlist = control_params$varlist)
    parallel::clusterExport(cl_use,
                            varlist = "grad_approx_fn",
                            envir = environment())
  }

  # pop initialization
  message('initalizing population...')
  if(control_params$parallel_type=='none'){
    # pop initialization sequentially
    for(pmem_index in 1:control_params$n_particles){
      count = 0 # establish a count variable to avoid infinite run time
      for(l in 1:length(ObjFun)){
        weights[1,pmem_index,l] <- Inf
        while(weights[1,pmem_index,l]==Inf) {

          # only sample the particles left to init
          particles[1, pmem_index, ] =
            stats::rnorm(control_params$n_params * length(pmem_index),
                         rep(x = control_params$init_center, each = length(pmem_index)),
                         rep(x = control_params$init_sd, each = length(pmem_index)))

          if(is.list(ObjFun)){
            # temp_weight <- 0
            for(l in seq_along(ObjFun)){
              weights[1,pmem_index,l] = ObjFun[[l]](particles[1, pmem_index, ], ...)
              # TODO: for blocking
              # like_weights[1,pmem_index,l] = weights[1,pmem_index,l] -
              # prior_function_list[[l]](particles[1, pmem_index, ], ...)
              if(!is.finite(weights[1, pmem_index, l])){
                weights[1, pmem_index, l] = Inf
              }
            }
          }else{
            stop("ObjFun is not iterable!")
          }
          count = count + 1
          if(count>control_params$give_up_init){
            stop('population initialization failed.
        inspect objective function or change init_center/init_sd to sample more
             likely parameter values')
          }
        }
        message(paste0(pmem_index, " / ", control_params$n_particles))
      }
    }
  } else {
    # pop initialization in parallel
    count = 0 # establish a count variable to avoid infinite run time
    weights[1,,] <- Inf # initialize to worst possible value
    while(any(!is.finite(weights[1,,]))) {
      # determine which indices still need init
      pmem_index <- which(weights[1,,1]==Inf)
      message(paste0(length(pmem_index), " / ", control_params$n_particles, " left"))
      # sample from init distribution
      particles[1, pmem_index, ] =
        stats::rnorm(control_params$n_params * length(pmem_index),
                     rep(x = control_params$init_center, each = length(pmem_index)),
                     rep(x = control_params$init_sd, each = length(pmem_index)))
      # parallel apply ObjFun on needed particles only
      for(l in seq_along(ObjFun)){
        weights[1, pmem_index, l] =
          parallel::parApply(cl = cl_use,
                             X = particles[1, pmem_index,] |>
                               matrix(nrow = control_params$n_particles),
                             MARGIN = c(1),
                             FUN = ObjFun[[l]], ...)
        # TODO: for blocking
        # like_weights[1,pmem_index,l] = weights[1,pmem_index,l] - parallel::parApply(cl = cl_use,
        # X = particles[1, pmem_index,],
        # MARGIN = c(1),
        # FUN = prior_function_list[[l]], ...)

        # catch NA's and Infinity and assign worst possible value
        if(any(!is.finite(weights[1, pmem_index, l]))){
          weights[1, any(!is.finite(weights[1, pmem_index, l])), ] = Inf
        }
      }
      count = count + 1
      if(count>control_params$give_up_init){
        stop('population initialization failed.
        inspect objective function or change init_center/init_sd to sample more
             likely parameter values')
      }
    }
  }
  message('population initialization complete  :)')

  message("running SQG-DE...")

  iter_idx=1
  converge_test_passed=FALSE
  for(iter in 1:control_params$n_iter){

    if(control_params$parallel_type=='none'){
      for(l in seq_along(ObjFun)){
        # adapt particles using SQG DE sequentially
        temp=matrix(unlist(lapply(1:control_params$n_particles,
                                  SQG_DE_bin_1_pos,
                                  current_params = particles[iter_idx,,] |>
                                    matrix(nrow = control_params$n_particles),   # current parameter values (numeric matrix)
                                  current_weight = weights[iter_idx,,l],  # corresponding weights (numeric vector)
                                  objFun = ObjFun[[l]],  # objective function (returns scalar)
                                  scheme = control_params$adapt_scheme, # TODO: later feature for conciseness
                                  # current_like_weight = like_weights[iter_idx,,l], # TODO: later update to only update the prior density in blocked updating
                                  # prior_function = prior_function_list[[l]],
                                  resample_weight = control_params$resample_weight,
                                  params_update_ind_vec = control_params$param_ind_to_update_list[[l]],
                                  step_size = control_params$step_size,
                                  jitter_size = control_params$jitter_size,
                                  n_particles = control_params$n_particles,
                                  crossover_rate = control_params$crossover_rate,
                                  n_diff = control_params$n_diff,
                                  ...)),
                    control_params$n_particles,
                    control_params$n_params+1, byrow=TRUE)

        # update particle after adaption
        weights[iter_idx,,l] = temp[, 1]
        # like_weights[iter_idx,,l] = temp[, 2] # TODO: later update to only update the prior density in blocked updating
        particles[iter_idx, , control_params$param_ind_to_update_list[[l]]] =
          temp[, c(FALSE, control_params$param_ind_to_update_list[[l]])]
        # particles[iter_idx, , control_params$param_ind_to_update_list[[l]]] =
        #   temp[, c(FALSE, FALSE, control_params$param_ind_to_update_list[[l]])]
      }
    } else {
      for(l in seq_along(ObjFun)){
        # adapt particles using SQG DE in parallel
        temp=matrix(unlist(parallel::parLapplyLB(cl = cl_use,
                                                 X = 1:control_params$n_particles,
                                                 fun = SQG_DE_bin_1_pos,
                                                 current_params = particles[iter_idx,,] |>
                                                   matrix(nrow = control_params$n_particles),   # current parameter values (numeric matrix)
                                                 current_weight = weights[iter_idx,,l],  # corresponding weights (numeric vector)
                                                 objFun = ObjFun[[l]],  # objective function (returns scalar)
                                                 scheme = control_params$adapt_scheme, # TODO: later feature for conciseness
                                                 # current_like_weight = like_weights[iter_idx,,l], # TODO: later update to only update the prior density in blocked updating
                                                 # prior_function = prior_function_list[[l]],
                                                 resample_weight = control_params$resample_weight,
                                                 params_update_ind_vec = control_params$param_ind_to_update_list[[l]],
                                                 step_size = control_params$step_size,
                                                 jitter_size = control_params$jitter_size,
                                                 n_particles = control_params$n_particles,
                                                 crossover_rate = control_params$crossover_rate,
                                                 n_diff = control_params$n_diff,
                                                 ...)),
                    control_params$n_particles,
                    control_params$n_params+1, byrow=TRUE)
        # update particle after adaption
        weights[iter_idx,,l] = temp[, 1]
        # like_weights[iter_idx,,l] = temp[, 2] # TODO: later update to only update the prior density in blocked updating
        particles[iter_idx, , control_params$param_ind_to_update_list[[l]]] =
          temp[, c(FALSE, control_params$param_ind_to_update_list[[l]])]
        # particles[iter_idx, , control_params$param_ind_to_update_list[[l]]] =
        #   temp[, c(FALSE, FALSE, control_params$param_ind_to_update_list[[l]])]
      }
    }
    # carry over particles for next iteration
    if(iter<control_params$n_iter){
      weights[iter_idx+1,,] = weights[iter_idx,,]
      particles[iter_idx+1,,] = particles[iter_idx,,]
    }

    #####################
    ####### purify
    #####################
    if(iter%%control_params$purify==0){

      if(control_params$parallel_type=='none'){
        for(l in seq_along(ObjFun)){
          temp=matrix(unlist(lapply(1:control_params$n_particles,
                                    Purify,
                                    current_params = particles[iter_idx,,] |>
                                      matrix(nrow = control_params$n_particles),   # current parameter values (numeric  matrix)
                                    params_update_ind_vec = control_params$param_ind_to_update_list[[l]],
                                    current_weight = weights[iter_idx,,l],  # corresponding weights (numeric vector)
                                    objFun = ObjFun[[l]],  # objective function (returns scalar)
                                    scheme = control_params$adapt_scheme, # TODO: later feature for conciseness
                                    # current_like_weight = like_weights[iter_idx,,l], # TODO: later update to only update the prior density in blocked updating
                                    # prior_function = prior_function_list[[l]],
                                    ...)),
                      nrow = control_params$n_particles,
                      ncol = control_params$n_params+1, byrow=TRUE)

          # update particle after adaption
          weights[iter_idx,,l] = temp[, 1]
          # like_weights[iter_idx,,l] = temp[, 2] # TODO: later update to only update the prior density in blocked updating
          particles[iter_idx, , control_params$param_ind_to_update_list[[l]]] =
            temp[, c(FALSE, control_params$param_ind_to_update_list[[l]])]
          # particles[iter_idx, , control_params$param_ind_to_update_list[[l]]] =
          #   temp[, c(FALSE, FALSE, control_params$param_ind_to_update_list[[l]])]
        }
      } else {
        for(l in seq_along(ObjFun)){
          temp=matrix(unlist(parallel::parLapplyLB(cl_use,
                                                   1:control_params$n_particles,
                                                   Purify,
                                                   current_params = particles[iter_idx,,] |>
                                                     matrix(nrow = control_params$n_particles),   # current parameter values (numeric matrix)
                                                   params_update_ind_vec = control_params$param_ind_to_update_list[[l]],
                                                   current_weight = weights[iter_idx,,l],  # corresponding weights (numeric vector)
                                                   objFun = ObjFun[[l]],  # objective function (returns scalar)
                                                   scheme = control_params$adapt_scheme, # TODO: later feature for conciseness
                                                   # current_like_weight = like_weights[iter_idx,,l], # TODO: later update to only update the prior density in blocked updating
                                                   # prior_function = prior_function_list[[l]],
                                                   ...)),
                      control_params$n_particles,
                      control_params$n_params+1, byrow=TRUE)
          # update particle after adaption
          weights[iter_idx,,l] = temp[, 1]
          # like_weights[iter_idx,,l] = temp[, 2] # TODO: later update to only update the prior density in blocked updating
          particles[iter_idx, , control_params$param_ind_to_update_list[[l]]] =
            temp[, c(FALSE, control_params$param_ind_to_update_list[[l]])]
          # particles[iter_idx, , control_params$param_ind_to_update_list[[l]]] =
          #   temp[, c(FALSE, FALSE, control_params$param_ind_to_update_list[[l]])]
        }
        # update particle after adaption
        weights[iter_idx, ] = temp[, 1]
        particles[iter_idx, , ] = temp[, 2:(control_params$n_params+1)]
        # carry over particles for next iteration
        if(iter<control_params$n_iter){
          weights[iter_idx+1, ] = temp[, 1]
          particles[iter_idx+1, , ] = temp[, 2:(control_params$n_params+1)]
        }
      }
    }

    if(iter %% control_params$stop_check==0){
      # assign convergence method
      if(control_params$converge_crit=='percent'){
        if(length(weights[1,1,]) > 1){
          percent_improve=(1-((weights[iter_idx,,, drop = FALSE] |> apply(MARGIN = 1, stats::median)))/
                             (weights[iter_idx-control_params$stop_check+1,,, drop = FALSE] |> apply(1, stats::median)))*100
          if(all(percent_improve<(control_params$stop_tol))){
            message("Convergence criterion met. Stopping optimization early")
            converge_test_passed=TRUE
            break
          }
        }else{
          percent_improve=(1-(weights[iter_idx,,1] |> stats::median())/
                             (weights[iter_idx-control_params$stop_check+1,,1] |> stats::median()))*100
          if(percent_improve<control_params$stop_tol){
            message("Convergence criterion met. Stopping optimization early")
            converge_test_passed=TRUE
            break
          }
        }
      }
      if(control_params$converge_crit=='stdev'){
        if(length(weights[1,1,]) > 1){
          weight_sd=(weights[iter_idx:(iter_idx-control_params$stop_check+1),,, drop = FALSE]) |>
            apply(1, stats::sd)
          if(all(weight_sd<(control_params$stop_tol))){
            message("Convergence criterion met. Stopping optimization early")
            converge_test_passed=TRUE
            break
          }
        }else{
          weight_sd=(weights[iter_idx:(iter_idx-control_params$stop_check+1),,1]) |>
            stats::sd()
          if(weight_sd<control_params$stop_tol){
            message("Convergence criterion met. Stopping optimization early")
            converge_test_passed=TRUE
            break
          }
        }
      }
    }


    SQGDE_out <<- list('particles_trace' = particles,
                       'weights_trace' = weights,
                       'converged' = converge_test_passed)


    if(iter%% control_params$iter_message_freq == 0){
      message(paste0('iter ', iter, '/', control_params$n_iter))
    }
    if(iter%%control_params$thin==0 & !(iter==control_params$n_iter)){
      iter_idx = iter_idx+1
    }
    if(control_params$save_int > 0 & iter%% control_params$save_int == 0){
      message(paste0('saving ouput on iter ', iter, '/', control_params$n_iter))
      saveRDS(SQGDE_out, file = control_params$save_rds_string)
      message(paste0('ouput saved'))
    }
  }
  message(paste0('run complete!'))
  # cluster stop
  if(!control_params$parallel_type=='none'){
    parallel::stopCluster(cl = cl_use)
  }
  minIdx = which.min(weights[iter_idx,,, drop = FALSE] |> apply(1, sum))
  minEst = particles[iter_idx, minIdx, ]

  if(control_params$return_trace==TRUE){
    return(list('solution' = minEst,
                'weight' = weights[iter_idx, minIdx,],
                'particles_trace' = particles,
                'weights_trace' = weights,
                'converged' = converge_test_passed))
  } else {
    return(list('solution' = minEst,
                'weight' = weights[iter_idx, minIdx,],
                'converged' = converge_test_passed))
  }
}
