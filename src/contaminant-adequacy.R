library(R2jags)
load(file = "data/memory.RData")
n_trials <- dim(results$res)[1]
n_participants <- dim(results$res)[2]
n_groups <- dim(results$res)[3]
n_st_type <- length(unique(results$st[,1,1]))

jags_data <- list(n_trials = n_trials,
                  n_participants = n_participants,
                  n_groups = n_groups,
                  stimulus = results$st,
                  responses = results$res)

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
  
# Posterior adequacy
  for(a in 1:n_groups){
    for(i in 1:n_participants){
      for(t in 1:n_trials){
        theta_hat[a,i,t] = theta[a,i,t,z[a,i,t]+1]
      }
    }
  }
}"

jags_parameter <- c("theta_hat")

samples <- jags(data = jags_data,
                parameters.to.save = jags_parameter,
                model.file = textConnection(jags_model), 
                n.chains = 4, 
                n.iter = 20000,
                n.burnin = 10000,
                n.thin = 1,
                DIC = TRUE)

theta_hat <- samples$BUGSoutput$mean$theta_hat

save(theta_hat,file='data/contaminant-mean-theta.RData')