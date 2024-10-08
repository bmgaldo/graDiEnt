% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetAlgoParams.R
\name{GetAlgoParams}
\alias{GetAlgoParams}
\title{GetAlgoParams}
\usage{
GetAlgoParams(
  n_params,
  n_particles = NULL,
  n_diff = 2,
  n_iter = 1000,
  init_sd = 0.01,
  init_center = 0,
  n_cores_use = 1,
  step_size = NULL,
  jitter_size = 1e-06,
  crossover_rate = 1,
  parallel_type = "none",
  return_trace = FALSE,
  thin = 1,
  purify = Inf,
  adapt_scheme = NULL,
  give_up_init = 100,
  stop_check = 10,
  stop_tol = 1e-04,
  converge_crit = "stdev",
  var_list = NULL,
  message_int = 100
)
}
\arguments{
\item{n_params}{The number of parameters estimated/optimized, this integer value NEEDS to be specified.}

\item{n_particles}{The number of particles (population size), 3*n_params is the default value.}

\item{n_diff}{The number of mutually exclusive vector pairs to stochastically approximate the gradient.}

\item{n_iter}{The number of iterations to run the algorithm, 1000 is default.}

\item{init_sd}{A positive scalar or n_params-dimensional numeric vector, determines the standard deviation of the Gaussian initialization distribution. The default value is 0.01.}

\item{init_center}{A scalar or n_params-dimensional numeric vector, determines the mean of the Gaussian initialization distribution. The default value is 0.}

\item{n_cores_use}{An integer specifying the number of cores used when using parallelization. The default value is 1.}

\item{step_size}{A positive scalar, jump size or "F" in the DE crossover step notation. The default value is 2.38/sqrt(2*n_params).}

\item{jitter_size}{A positive scalar that determines the jitter (noise) size. Noise is added during adaption step from Uniform(-jitter_size,jitter_size) distribution. 1e-6 is the default value. Set to 0 to turn off jitter.}

\item{crossover_rate}{A numeric scalar on the interval (0,1]. Determines the probability a parameter on a chain is updated on a given crossover step, sampled from a Bernoulli distribution. The default value is 1.}

\item{parallel_type}{A character string specifying parallelization type. 'none','FORK', or 'PSOCK' are valid values. 'none' is default value. 'FORK' does not work with Windows OS.}

\item{return_trace}{A boolean, if true, the function returns particle trajectories. This is helpful for assessing convergence or debugging model code. The trace will be an iteration/thin $x$ n_particles $x$ n_params array containing parameter values and an iteration/thin $x$ n_particles array containing particle weights.}

\item{thin}{A positive integer. Only every 'thin'-th iteration will be stored in memory. The default value is 1. Increasing thin will reduce the memory required when running the algorithim for longer.}

\item{purify}{A positive integer. On every 'purify'-th iteration the particle weights are recomputed. This is useful if the objective function is stochastic/noisy. If the objective function is deterministic, this computation is redundant. Purify is set to Inf by default, disabling it.}

\item{adapt_scheme}{A character string that must be 'rand','current', or 'best' that determines the DE adaption scheme/strategy. 'rand' uses rand/1/bin DE-like scheme where a random particle and the particle-based quasi-gradient approximation are used to generate proposal updates for a given particle. 'current' uses current/1/bin, and 'best' uses best/1/bin which follow an analogous adaption scheme to rand. 'rand' is the default value.}

\item{give_up_init}{An integer for how many failed initialization attempts before stopping the optimization routine. 100 is the default value.}

\item{stop_check}{An integer for how often to check the convergence criterion. The default is 10 iterations.}

\item{stop_tol}{A convergence metric must be less than value to be labeled as converged. The default is 1e-4.}

\item{converge_crit}{A character string denoting the convergence metric used, valid metrics are 'stdev' (standard deviation of population weight in the last stop_check iterations) and 'percent' (percent improvement in median particle weight in the last stop_check iterations). 'stdev' is the default.}

\item{var_list}{A vector of names of variables and functions to export for parallelization (via parallel::clusterExport). Default value is NULL.}

\item{message_int}{An integer specifying the interval at which an update message is displayed, based on the number of completed iterations. Must be a positive value. The default is 100; set to \code{Inf} to disable these particular messages.}
}
\value{
A list of control parameters for the optim_SQGDE function.
}
\description{
Get control parameters for optim_SQGDE function.
}
