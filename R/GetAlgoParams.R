#' GetAlgoParams
#'
#' @description Generates a list of control parameters for the `optim_SQGDE` function, with defaults and validation checks for robustness and ease of use.
#'
#' @param n_params **(Required)** Integer. The number of parameters to be optimized. Must be a positive scalar.
#' @param param_ind_to_update_list List of vectors indicating which parameter indices to update for each sub-function in `ObjFun`. For example, `param_ind_to_update_list[[L]]` specifies the subset of parameters updated for the `L`-th function. Default is `NULL`, which updates all parameters for all functions.
#' @param resample_weight Logical. If `TRUE`, resamples weights for parameters based on changes since the last weight evaluation. Default is `FALSE`.
#' @param n_particles Integer. Number of particles (population size). Default is `3 * n_params`.
#' @param n_iter Integer. Number of iterations for the algorithm. Default is `1000`.
#' @param n_diff Integer. Number of mutually exclusive vector pairs used to approximate the gradient. Default is `2`.
#' @param crossover_rate Numeric in (0, 1]. Probability of updating a parameter during a crossover step. Default is `1`.
#' @param init_sd Numeric. Standard deviation for the Gaussian initialization distribution. Can be a scalar or vector of length `n_params`. Default is `0.01`.
#' @param init_center Numeric. Mean for the Gaussian initialization distribution. Can be a scalar or vector of length `n_params`. Default is `0`.
#' @param n_cores_use Integer. Number of cores for parallelization. Default is `1`.
#' @param step_size Numeric. Jump size for the differential evolution crossover step. Default is `2.38 / sqrt(2 * n_params)`.
#' @param jitter_size Numeric. Adds uniform noise (`Uniform(-jitter_size, jitter_size)`) during the adaptation step. Default is `1e-6`. Set to `0` to disable.
#' @param parallel_type Character. Type of parallelization: `'none'`, `'FORK'`, or `'PSOCK'`. Default is `'none'`. Note: `'FORK'` does not work on Windows.
#' @param parallel_seed Integer. Seed for reproducibility in parallel clusters. Default is `NULL`.
#' @param return_trace Logical. If `TRUE`, the function stores and returns particle trajectories for debugging or convergence assessment. Default is `FALSE`.
#' @param thin Integer. Thinning interval for storing iterations. Default is `1`. Higher values reduce memory usage.
#' @param purify Integer or `Inf`. Interval for recomputing particle weights. Useful for stochastic objective functions. Default is `Inf` (disabled).
#' @param adapt_scheme Character. Adaptation scheme: `'rand'`, `'current'`, or `'best'`. Default is `'rand'`.
#' @param give_up_init Integer. Number of failed initialization attempts before stopping. Default is `100`.
#' @param stop_check Integer. Frequency (in iterations) of convergence checks. Default is `10`.
#' @param stop_tol Numeric. Convergence tolerance for stopping criteria. Default is `1e-4`.
#' @param converge_crit Character. Convergence metric: `'stdev'` (standard deviation of weights) or `'percent'` (percent improvement in median weight). Default is `'stdev'`.
#' @param varlist List. Variables and functions to export during parallelization.
#' @param iter_message_freq Integer. Frequency (in iterations) of console messages. Default is `10`.
#' @param save_int Integer. Interval (in iterations) for saving progress. Negative values disable saving. Default is `-1`.
#' @param save_rds_string Character. File name for saving function output if `save_int > 0`. Default is `"SQGDE_DEFAULT_IMAGE.rds"`.
#' @param outfile_string Character. File name for console output logging. Default is `"console_output.txt"`.
#'
#' @return A list of validated control parameters for the `optim_SQGDE` function.
#' @export
GetAlgoParams = function(n_params,
                         param_ind_to_update_list = NULL,
                         resample_weight = FALSE,
                         n_particles = NULL,
                         n_diff = 2,
                         n_iter = 1000,
                         init_sd = 0.01,
                         init_center = 0,
                         n_cores_use = 1,
                         step_size = NULL,
                         jitter_size = 1e-6,
                         crossover_rate = 1,
                         parallel_type = 'none',
                         parallel_seed = NULL,
                         return_trace = FALSE,
                         thin = 1,
                         purify = NULL,
                         adapt_scheme = NULL,
                         give_up_init = 100,
                         stop_check = 10,
                         stop_tol = 1e-4,
                         converge_crit = 'stdev',
                         varlist = NULL,
                         iter_message_freq = 100,
                         save_int = -1,
                         save_rds_string = "SQGDE_DEFAULT_IMAGE.rds"){
  # n_params
  ### catch errors
  if(length(n_params) != 1){
    stop('ERROR: n_params is not a finite number')
  }else if (!is.numeric(n_params) | !is.finite(n_params)){
    stop('ERROR: n_params is not a finite number')
  }else{
    n_params = as.integer(n_params)
  }
  if(any(!is.finite(n_params))){
    stop('ERROR: n_params is not finite')
  }  else if( n_params<1 | length(n_params)>1){
    stop('ERROR: n_params must be a postitive integer scalar')
  }

  ### no list given, default to all parameters
  if(is.null(param_ind_to_update_list)){
    param_ind_to_update_list = list(rep(TRUE, n_params))
  }
  if(is.list(param_ind_to_update_list)){
    for(l in param_ind_to_update_list){
      if(!is.numeric(l)){
        if(!is.logical(l)){
          stop('ERROR: each element of param_ind_to_update_list must be a LOGICAL
             vector of parameter indices of length n_params')
        }
      }
      # if(is.numeric(l)){
      #   if(any(l<1) | any(length(l)>n_params)){
      #     stop('ERROR: value in the numeric vector of each element
      #        param_ind_to_update_list of must be a between 1 and n_params')
      #   }
      # }
      if(length(l) != n_params){
        stop('ERROR: each element of param_ind_to_update_list must be length
           n_params')
      }
      if(any(!is.finite(l))){
        stop('ERROR: each element param_ind_to_update_list of must be a finite
           numeric vector')
      }
    }
  }else{
    if(length(param_ind_to_update_list) != n_params){
      stop('ERROR: param_ind_to_update_list must be length
           n_params if not a list')
    }
    if(any(!is.finite(param_ind_to_update_list))){
      stop('ERROR: param_ind_to_update_list of must be a finite
           numeric vector if not a list')
    }
  }
  # if(length(param_ind_to_update_list) != length(ObjFun)){
  #   stop('ERROR: param_ind_to_update_list must be the same length as ObjFun')
  # }

  # resample_weight
  if(is.null(resample_weight)){
    resample_weight = FALSE
  }else{
    resample_weight = as.logical(resample_weight)
  }

  # n_particles
  ### if null assign default value
  if(is.null(n_particles)){
    n_particles = max(3*n_params,4)
  }
  ### catch errors
  n_particles = as.integer(n_particles)
  if(any(!is.finite(n_particles))){
    stop('ERROR: n_particles is not finite')
  } else if( n_particles<4 | length(n_particles)>1){
    stop('ERROR: n_particles must be a postitive integer scalar, and atleast 4')
  }

  # n_iter
  ### if null assign default value
  if(is.null(n_iter)){
    n_iter = 1000
  }
  ### catch errors
  n_iter = as.integer(n_iter)
  if(any(!is.finite(n_iter))){
    stop('ERROR: n_iter is not finite')
  } else if( n_iter<4 | length(n_iter)>1){
    stop('ERROR: n_iter must be a postitive integer scalar, and atleast 4')
  }

  # init_sd
  init_sd = as.numeric(init_sd)
  if(any(!is.finite(init_sd))){
    stop('ERROR: init_sd is not finite')
  } else if(any(init_sd<= 0 | is.complex(init_sd))){
    stop('ERROR: init_sd must be positive and real-valued')
  } else if(!(length(init_sd) == 1 | length(init_sd) == n_params)){
    stop('ERROR: init_sd vector length must be 1 or n_params')
  }

  # init_center
  init_center = as.numeric(init_center)
  if(any(!is.finite(init_center))){
    stop('ERROR: init_center is not finite')
  } else if(any(is.complex(init_center))){
    stop('ERROR: init_center must be real valued')
  } else if(!(length(init_center) == 1 | length(init_center) == n_params)){
    stop('ERROR: init_center vector length must be 1 or n_params')
  }

  # n_cores_use
  ### assign NULL value default
  if(is.null(n_cores_use)){
    n_cores_use = 1
  }
  ### catch any errors
  n_cores_use = as.integer(n_cores_use)
  if(any(!is.finite(n_cores_use))){
    stop('ERROR: n_cores_use is not finite')
  } else if( n_cores_use<1 | length(n_cores_use)>1){
    stop('ERROR: n_cores_use must be a postitive integer scalar, and atleast 1')
  }


  # step_size
  ### assign NULL value default
  if(is.null(step_size)){
    step_size = 2.38/sqrt(2*n_params) # step size recommend in ter braak's 2006 paper
  }
  ### catch any errors
  if(any(!is.finite(step_size))){
    stop('ERROR: step_size is not finite')
  } else if(any(step_size<= 0 | is.complex(step_size))){
    stop('ERROR: step_size must be positive and real-valued')
  } else if(!(length(step_size) == 1)){
    stop('ERROR: step_size vector length must be 1 ')
  }

  #jitter_size
  ### assign NULL value default
  if(is.null(jitter_size)){
    jitter_size = 1e-6
  }
  ### catch any errors
  if(any(!is.finite(jitter_size))){
    stop('ERROR: jitter_size is not finite')
  } else if(any(jitter_size<= 0 | is.complex(jitter_size))){
    stop('ERROR: jitter_size must be positive and real-valued')
  } else if(!(length(jitter_size) == 1)){
    stop('ERROR: jitter_size vector length must be 1 ')
  }

  # crossover_rate
  ### if null assign default value
  if(any(is.null(crossover_rate))){
    crossover_rate = 1
  }
  ### catch errors
  crossover_rate = as.numeric(crossover_rate)
  if(any(!is.finite(crossover_rate))){
    stop('ERROR: crossover_rate is not finite')
  } else if(any(crossover_rate>1) | any(crossover_rate<= 0) | length(crossover_rate)>1){
    stop('ERROR: crossover_rate must be a numeric scalar on the interval (0,1]')
  } else if(is.complex(crossover_rate)){
    stop('ERROR: crossover_rate cannot be complex')
  }

  #parallel_type
  validParType = c('none','FORK','PSOCK')
  ### assign NULL value default
  if(is.null(parallel_type)){
    parallel_type = 'none'
  }
  ### catch any errors
  if(!parallel_type %in% validParType){
    stop(paste('ERROR: invalid parallel_type.'))
  }

  #parallel_seed
  ### assign NULL value default
  if(is.null(parallel_seed)){
    parallel_seed = NULL
  }else if(!(parallel_seed %% 1 == 0) | parallel_seed < 0){
    ### catch any errors
    stop(paste('ERROR: invalid parallel_seed.
               Please select an integer greater than 0.'))
  }

  #converge_crit
  validConType = c('stdev','percent')
  ### assign NULL value default
  if(is.null(converge_crit)){
    converge_crit = 'stdev'
  }
  ### catch any errors
  if(!converge_crit %in% validConType){
    stop(paste('ERROR: invalid converge_crit.'))
  }

  validAdaptType = c('rand','current','best')
  ### assign NULL value default
  if(is.null(adapt_scheme)){
    adapt_scheme = 'rand'
  }
  ### catch any errors
  if(!adapt_scheme %in% validAdaptType){
    stop(paste('ERROR: invalid adaption scheme.'))
  }

  # thin
  ### if null assign default value
  if(is.null(thin)){
    thin = 1
  }
  ### catch errors
  thin = as.integer(thin)
  if(any(!is.finite(thin))){
    stop('ERROR: thin is not finite')
  } else if(any(thin<1) | length(thin)>1){
    stop('ERROR: thin must be a scalar postive integer')
  }

  #number iters stored per particle
  n_iters_per_particle = floor((n_iter)/thin)
  ### catch errors
  if(n_iters_per_particle<1 | (!is.finite(n_iters_per_particle))){
    stop('ERROR: number of stored particle value is negative or non finite.
         n_iters_per_particle = floor((n_iter)/thin)')
  }

  # purify
  no_check <- FALSE
  if(is.null(purify)){
    purify = Inf
    no_check <- TRUE
  }else if(!is.finite(purify) & !no_check){
    warning('Warning: purify was not finite. Defaulting to not using Purify.')
    purify = Inf
  }else if(is.finite(purify)){
    purify = as.integer(purify)
  }
  ### catch errors
  if(any(purify < 1) | (length(purify)>1)){
    stop('ERROR: purify must be a positive integer or Inf')
  }

  ### catch errors
  n_diff = as.integer(n_diff)
  if(any(!is.finite(n_diff))){
    stop('ERROR: n_diff is not finite')
  } else if(any(n_diff<1) | length(n_diff)>1){
    stop('ERROR: n_diff must be a scalar postive integer')
  } else if(n_diff>(n_particles/2)){
    stop('ERROR: n_diff cannot exceed n_particles/2')
  }

  ##################
  # give_up_init
  if(any(is.null(give_up_init)) | any(!is.finite(give_up_init))){
    warning('Warning: give_up_init was null or not finite. Using default value of 100.')
    give_up_init = 100
  }
  give_up_int=round(give_up_init)
  ### catch errors
  if(any(give_up_init < 1) | (length(give_up_init)>1)){
    stop('ERROR: give_up_init must be a scalar positive integer')
  }

  ##################
  # stop_check
  if(any(is.null(stop_check)) | any(!is.finite(stop_check))){
    warning('Warning: stop_check was null or non-finite. Using default value of 10.')
    stop_check = 10
  }
  stop_check=round(stop_check)
  ### catch errors
  if(any(stop_check < 2) | (length(stop_check)>1)){
    stop('ERROR: stop_check must be a scalar positive integer and greater than 1')
  }


  ##################
  # stop_tol
  if(any(is.null(stop_tol)) | any(is.na(stop_tol))){
    stop_tol = 1e-4
  }
  ### catch errors
  if(any(stop_tol < 0) | (length(stop_tol)>1)){
    stop('ERROR: stop_tol must be a scalar positive')
  }

  ##################
  # varlist
  if(any(!is.list(varlist))){
    varlist = as.list(varlist)
  }
  ### catch errors
  if(!is.list(varlist)){
    stop('ERROR: varlist must be a list of strings of variable and function name!')
  }
  for(l in varlist){
    if(!is.character(l)){
      stop('ERROR: each element of stop varlist must be a string of variable and function name!')
    }
  }

  ##################
  # iter_message_freq
  if(is.null(iter_message_freq)){
    iter_message_freq = 100
  }
  ### catch errors
  iter_message_freq = as.integer(iter_message_freq)
  if(any(!is.finite(iter_message_freq))){
    stop('ERROR: iter_message_freq is not finite')
  } else if( iter_message_freq<1 | length(iter_message_freq)>1){
    stop('ERROR: iter_message_freq must be a postitive integer scalar, and atleast 1')
  }

  ##################
  # save_int
  no_save_flag <- FALSE
  if(is.null(save_int)){
    save_int = -1
    no_save_flag <- TRUE
  }else if(!is.finite(save_int)){
    stop('ERROR: save_int is not finite')
  }else if(save_int < 0){
    save_int = -1
    no_save_flag <- TRUE
  }
  ### catch errors
  save_int = as.integer(save_int)
  if(any(!is.finite(save_int))){
    stop('ERROR: save_int is not finite')
  }
  if(!no_save_flag){
    if(length(save_int) > 1 | save_int < 1){
      stop('ERROR: save_int must be a postitive integer scalar, and at least 1 or
         less than 0 to not save')
    }
  }



  ##################
  # save_rds_string
  if(save_int > 0){
    if(!is.character(save_rds_string)){
      stop('ERROR: save_rds_string must be a single CHARACTER sting that ends in .rds')
    }
    if(length(save_rds_string) != 1){
      stop('ERROR: save_rds_string must be a SINGLE character sting that ends in .rds')
    }
    if(substring(text = save_rds_string,
                 nchar(save_rds_string)-3,
                 nchar(save_rds_string)) != ".rds"){
      warning('Warning: save_rds_string must be a character sting that ends in .rds.
              Adding .rds to supplied string.')
      save_rds_string <- paste0(save_rds_string, ".rds")
    }
  }


  ##################
  # outfile_string
  if(!is.character(outfile_string)){
      stop('ERROR: outfile_string must be a single CHARACTER sting that ends in .txt')
  }
  if(length(outfile_string) != 1){
      stop('ERROR: outfile_string must be a SINGLE character sting that ends in .txt')
  }
  if(substring(text = outfile_string,
                 nchar(outfile_string)-3,
                 nchar(outfile_string)) != ".txt"){
      warning('Warning: outfile_string must be a character sting that ends in .txt.
              Adding .txt to supplied string.')
      outfile_string <- paste0(outfile_string, ".txt")
  }

  out = list('n_params' = n_params,
             'param_ind_to_update_list' = param_ind_to_update_list,
             'resample_weight' = resample_weight,
             'n_particles' = n_particles,
             'n_iter' = n_iter,
             'init_sd' = init_sd,
             'init_center' = init_center,
             'n_cores_use' = n_cores_use,
             'step_size' = step_size,
             'crossover_rate' = crossover_rate,
             'jitter_size' = jitter_size,
             'parallel_type' = parallel_type,
             'thin' = thin,
             'purify' = purify,
             'n_iters_per_particle' = n_iters_per_particle,
             'return_trace' = return_trace,
             'n_diff' = n_diff,
             'adapt_scheme' = adapt_scheme,
             'give_up_init'= give_up_init,
             'stop_tol' = stop_tol,
             'stop_check' = stop_check,
             'converge_crit' = converge_crit,
             'varlist' = varlist,
             'iter_message_freq' = iter_message_freq,
             'parallel_seed' = parallel_seed,
             'save_int' = save_int,
             'save_rds_string' = save_rds_string,
             'outfile' = outfile_string)

  return(out)
}
