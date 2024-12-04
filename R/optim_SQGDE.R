#' optim_SQGDE
#'
#' @description  Runs Stochastic Quasi-Gradient Differential Evolution (SQG-DE;
#' Sala, Baldanzini, and Pierini, 2018) to minimize an objective function f(x).
#' To maximize a function f(x), simply pass g(x)=-f(x) to ObjFun argument.
#' @param ObjFun_list A list or sinngle function that is a scalar-returning
#' function(s) to minimize whose first argument is a real-valued
#' n_params-dimensional vector. A list of functions enables "blocking", where
#' the list of functions are "sub-iterated" through so that only a subset of
#' parameters are updated at one sub-iteration. Currently, this feature requires
#' the user to specify the "Master Objective Function" where all particles are
#' updated as the first element of the list for initialization.
#' @param control_params control parameters for SQG-DE algo. see
#' \code{\link{GetAlgoParams}} function documentation for more details. The only
#' argument you NEED to pass here is n_params.
#' @param ... additional arguments to pass ObjFun.
#' @return list containing solution and it's corresponding weight (i.e. f(solution)).
#' @export
#' @md
#' @examples
#' ##############
#' # Maximum Likelihood Example
#' ##############
#'
#' # simulate from model
#' dataExample=matrix(rnorm(1000,c(-1,1,0,1),c(1,1,1,1)),ncol=4,byrow = TRUE)
#'
#' # list parameter names
#' param_names_example=c("mu_1","mu_2","mu_3","mu_4")
#'
#' # negative log likelihood
#' ExampleObjFun=function(x,data,param_names){
#'   out=0
#'
#'   names(x) <- param_names
#'
#'   # log likelihoods
#'   out=out+sum(dnorm(data[,1],x["mu_1"],sd=1,log=TRUE))
#'   out=out+sum(dnorm(data[,2],x["mu_2"],sd=1,log=TRUE))
#'   out=out+sum(dnorm(data[,3],x["mu_3"],sd=1,log=TRUE))
#'   out=out+sum(dnorm(data[,4],x["mu_4"],sd=1,log=TRUE))
#'
#'   return(out*-1)
#' }
#'
#' ########################
#' # run optimization
#' out <- optim_SQGDE(ObjFun_list = ExampleObjFun,
#'                    control_params = GetAlgoParams(n_params = length(param_names_example),
#'                                              n_iter = 250,
#'                                               n_particles = 12,
#'                                               n_diff = 2,
#'                                               return_trace = TRUE),
#'                    data = dataExample,
#'                    param_names = param_names_example)
#' old_par <- par() # save graphic state for user
#' # plot particle trajectory
#'
#' par(mfrow=c(2,2))
#' matplot(out$particles_trace[,,1],type='l')
#' matplot(out$particles_trace[,,2],type='l')
#' matplot(out$particles_trace[,,3],type='l')
#' matplot(out$particles_trace[,,4],type='l')
#'
#' #SQG DE solution
#' out$solution
#'
#' #analytic solution
#' apply(dataExample, 2, mean)
#'
#' par(old_par) # restore user graphic state
#'

optim_SQGDE = function(ObjFun_list,
                       # prior_function_list = NULL,
                       control_params = GetAlgoParams(), ...){


  # ObjFun_list
  ### catch errors
  if(length(ObjFun_list) == 1){
    if(!is.function(ObjFun_list)){
      stop('ERROR: ObjFun_list is not A FUNCTION or list of functions!')
    }
    ObjFun_list <- list(ObjFun_list)
  }else{
    if(!is.list(ObjFun_list)){
      stop('ERROR: ObjFun_list is not a function or LIST of functions!')
    }
    for(f in seq_along(ObjFun_list)){
      if(!is.function(ObjFun_list[[f]])){
        stop(paste0('ERROR: The ', f, 'element of ObjFun_list is not a function!'))
      }
    }
  }
  if(length(ObjFun_list) != length(control_params$param_ind_to_update_list)){
    stop(paste0('ERROR: The element of ObjFun_list has a length of ',
                length(ObjFun_list),
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
  #   if(length(prior_function_list) != length(ObjFun_list)){
  #     stop(paste0('ERROR: The prior_function_list is not the same length as ObjFun_list!'))
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
                          ObjFun_list |> length()))
  # TODO: for blocking
  # like_weights = array(NA,
  #                      dim = c(control_params$n_iters_per_particle,
  #                              control_params$n_particles,
  #                              ObjFun_list |> length()))

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
          outfile = "SQGDE_parallel_outfile.txt")
    }else if(control_params$parallel_type == "FORK"){
      cl_use =
        parallel::makeForkCluster(control_params$n_cores_use,
                                  outfile = "SQGDE_parallel_outfile.txt")
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
      for(l in 1:length(ObjFun_list)){
        weights[1,pmem_index,l] <- Inf
        while(weights[1,pmem_index,l]==Inf) {

          # only sample the particles left to init
          particles[1, pmem_index, ] =
            stats::rnorm(control_params$n_params * length(pmem_index),
                         rep(x = control_params$init_center, each = length(pmem_index)),
                         rep(x = control_params$init_sd, each = length(pmem_index)))

          if(is.list(ObjFun_list)){
            # temp_weight <- 0
            for(l in seq_along(ObjFun_list)){
              weights[1,pmem_index,l] = ObjFun_list[[l]](particles[1, pmem_index, ], ...)
              # TODO: for blocking
              # like_weights[1,pmem_index,l] = weights[1,pmem_index,l] -
              # prior_function_list[[l]](particles[1, pmem_index, ], ...)
              if(!is.finite(weights[1, pmem_index, l])){
                weights[1, pmem_index, l] = Inf
              }
            }
          }else{
            stop("ObjFun_list is not iterable!")
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
      for(l in seq_along(ObjFun_list)){
        weights[1, pmem_index, l] =
          parallel::parApply(cl = cl_use,
                             X = particles[1, pmem_index,] |>
                               matrix(nrow = control_params$n_particles),
                             MARGIN = c(1),
                             FUN = ObjFun_list[[l]], ...)
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
      for(l in seq_along(ObjFun_list)){
        # adapt particles using SQG DE sequentially
        temp=matrix(unlist(lapply(1:control_params$n_particles,
                                  SQG_DE_bin_1_pos,
                                  current_params = particles[iter_idx,,] |>
                                    matrix(nrow = control_params$n_particles),   # current parameter values (numeric matrix)
                                  current_weight = weights[iter_idx,,l],  # corresponding weights (numeric vector)
                                  objFun = ObjFun_list[[l]],  # objective function (returns scalar)
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
      for(l in seq_along(ObjFun_list)){
        # adapt particles using SQG DE in parallel
        temp=matrix(unlist(parallel::parLapplyLB(cl = cl_use,
                                                 X = 1:control_params$n_particles,
                                                 fun = SQG_DE_bin_1_pos,
                                                 current_params = particles[iter_idx,,] |>
                                                   matrix(nrow = control_params$n_particles),   # current parameter values (numeric matrix)
                                                 current_weight = weights[iter_idx,,l],  # corresponding weights (numeric vector)
                                                 objFun = ObjFun_list[[l]],  # objective function (returns scalar)
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
        for(l in seq_along(ObjFun_list)){
          temp=matrix(unlist(lapply(1:control_params$n_particles,
                                    Purify,
                                    current_params = particles[iter_idx,,] |>
                                      matrix(nrow = control_params$n_particles),   # current parameter values (numeric  matrix)
                                    params_update_ind_vec = control_params$param_ind_to_update_list[[l]],
                                    current_weight = weights[iter_idx,,l],  # corresponding weights (numeric vector)
                                    objFun = ObjFun_list[[l]],  # objective function (returns scalar)
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
        for(l in seq_along(ObjFun_list)){
          temp=matrix(unlist(parallel::parLapplyLB(cl_use,
                                                   1:control_params$n_particles,
                                                   Purify,
                                                   current_params = particles[iter_idx,,] |>
                                                     matrix(nrow = control_params$n_particles),   # current parameter values (numeric matrix)
                                                   params_update_ind_vec = control_params$param_ind_to_update_list[[l]],
                                                   current_weight = weights[iter_idx,,l],  # corresponding weights (numeric vector)
                                                   objFun = ObjFun_list[[l]],  # objective function (returns scalar)
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


    if(iter%% control_params$print_int == 0){
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
