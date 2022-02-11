# Data analysis using the contaminant model in other to draw samples from the 
# posterior predictive distribution of P(R = 'new'). The posterior distribution
# is used to calculate the accuracy of the model using the posterior mode

# call R2jags library
library(R2jags)

# load experimental data
load(file = "data/memory.RData")

# number of trials per subject
n_trials <- dim(results$res)[1]

# number of participants per group
n_participants <- dim(results$res)[2]

# number of groups
n_groups <- dim(results$res)[3]

# number of stimulus types
n_st_type <- length(unique(results$st[,1,1]))

# data to be parsed to JAGS
jags_data <- list(n_trials = n_trials,
                  n_participants = n_participants,
                  n_groups = n_groups,
                  stimulus = results$st,
                  responses = results$res)

# save contaminant model as text object in R
jags_model <-"model{
# Prior for discriminability individual level
  for(a in 1:n_groups){
    for(i in 1:n_participants){
      for(j in 1:6){
        dtmp[a,i,j] ~ dnorm(0,1/3^2)T(0,)
      }
      d[a,i,1] = 0
      d[a,i,2:7] = sort(dtmp[a,i,])
    }
  }

# Prior hyperparameters and pchoice criterion
  for(a in 1:n_groups){
    for(i in 1:n_participants){
      k[a,i] ~ dnorm(0,1/4)
      phi[a,i] ~ dbeta(1,10) 
    }
  }
  
# Likelihood
  for(a in 1:n_groups){
    for(i in 1:n_participants){
      for(t in 1:n_trials){
        z[a,i,t] ~ dbern(phi[a,i])
        theta[a,i,t,1] = 1-phi(k[a,i]-d[a,i,stimulus[t,i,a]])
        theta[a,i,t,2] = 0.5
        responses[t,i,a] ~ dbern(theta[a,i,t,z[a,i,t]+1])
      }
    }
  }
  
# Posterior distribution of P(R = 'new')
  for(a in 1:n_groups){
    for(i in 1:n_participants){
      for(t in 1:n_trials){
        theta_hat[a,i,t] = theta[a,i,t,z[a,i,t]+1]
      }
    }
  }
}"

# posterior distributions to be saved in samples
jags_parameter <- c("theta_hat")

# save posterior distributions of contaminant model in samples object
samples <- jags(data = jags_data,
                parameters.to.save = jags_parameter,
                model.file = textConnection(jags_model), 
                n.chains = 4, 
                n.iter = 20000,
                n.burnin = 10000,
                n.thin = 1,
                DIC = TRUE)

# extract posterior distribution of P(R = 'new')
theta_hat <- samples$BUGSoutput$mean$theta_hat

# save posterior distributions of P(R = 'new')
save(theta_hat,file='data/contaminant-mean-theta.RData')