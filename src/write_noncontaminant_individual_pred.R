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

  # Likelihood
    for(t in 1:trials){
      theta[t] = 1-phi(k-d[stimulus[t]])
      y[t] ~ dbern(theta[t])
    }
  }
  ","models/noncontaminant_individual_pred.txt")