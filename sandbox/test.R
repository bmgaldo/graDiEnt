# simulate from model
dataExample=matrix(rnorm(1000,c(-1,1,0,1),c(.1,.1,.1,.1)),ncol=4,byrow = TRUE)

# list parameter names
param_names_example=c("mu_1","mu_2","mu_3","mu_4")

# negative log likelihood
ExampleObjFun=function(x,data,param_names){
  out=0

  names(x) <- param_names

  # log likelihoods
  out=out+sum(dnorm(data[,1],x["mu_1"],sd=.1,log=TRUE))
  out=out+sum(dnorm(data[,2],x["mu_2"],sd=.1,log=TRUE))
  out=out+sum(dnorm(data[,3],x["mu_3"],sd=.1,log=TRUE))
  out=out+sum(dnorm(data[,4],x["mu_4"],sd=.1,log=TRUE))

  return(out*-1)
}

########################
# run optimization
out <- optim_SQGDE(ObjFun = ExampleObjFun,
                   control_params = GetAlgoParams(n_params=length(param_names_example),
                                             n_iter=200,
                                              n_particles=12,
                                              n_diff=2,
                                              return_trace = TRUE),
                   data = dataExample,
                   param_names = param_names_example)

# plot particle trajectory
par(mfrow=c(2,2))
matplot(out$particles_trace[,,1],type='l')
matplot(out$particles_trace[,,2],type='l')
matplot(out$particles_trace[,,3],type='l')
matplot(out$particles_trace[,,4],type='l')

#SQG DE solution
out$solution

#analytic solution
apply(dataExample, 2, mean)
