# This function writes the contaminant model used when the data is not expanded 
# and is used at an individual level
write("
  model{
# Prior for discriminability individual level
  for(j in 1:6){
    dtmp[j] ~ dnorm(0,1/3^2)T(0,)
  }
  d[1] = 0
  d[2:7] = sort(dtmp[])

# Prior hyperparameters and pchoice criterion
  k ~ dnorm(0,1/4)
  phi ~ dbeta(1,10) 
  
# Likelihood
  for(t in 1:trials){
    z[t] ~ dbern(phi)
    theta[t,1] = 1-phi(k-d[stimulus[t]])
    theta[t,2] = 0.5
    y[t] ~ dbern(theta[t,z[t]+1])
  }
  
# Posterior adequacy
  for(t in 1:trials){
    theta_hat[t] = theta[t,z[t]+1]
  }
}
","models/contaminant_individual_pred.txt")
