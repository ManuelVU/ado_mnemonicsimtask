# Model comparison between contaminant and non-contaminant model by SD on 
# rate of contaminant behavior.

library(R2jags)

load(file = "data/memory.RData")
n_trials <- dim(results$res)[1]
n_participants <- dim(results$res)[2]
n_groups <- dim(results$res)[3]
n_st_type <- length(unique(results$st[,1,1]))

jags_data <- list(n_trials = n_trials,
                  n_participants = n_participants,
                  n_groups = n_groups,
                  n_st_type = n_st_type,
                  stimulus = results$st,
                  responses = results$res)

jags_model <- "model{
# Prior hyperparameters for discriminability
  for(g in 1:n_groups){
    for(i in 1:n_participants){
      for(j in 1:(n_st_type-1)){
        dtmp[j,i,g] ~ dnorm(0,1/3^2)T(0,)
      }
      d[i,1,g] = 0
      d[i,2:7,g] = sort(dtmp[,i,g])
    }
  }

# Prior hyperparameters and pchoice criterion
  for(g in 1:n_groups){
    for(i in 1:n_participants){
      k[i,g] ~ dnorm(0,1/4)
      phi[i,g] ~ dbeta(1,10) 
    }
  }

  # Likelihood
  for(g in 1:n_groups){
    for(i in 1:n_participants){
      for(t in 1:n_trials){
        z[t,i,g] ~ dbern(phi[i,g])
        theta[t,i,g,1] = 1-phi(k[i,g]-d[i,stimulus[t,i,g],g])
        theta[t,i,g,2] = 0.5
        responses[t,i,g] ~ dbern(theta[t,i,g,z[t,i,g]+1])
      }
    }
  }
}
"

jags_parameter <- c("d", "phi")

samples <- jags(data = jags_data,
                parameters.to.save = jags_parameter,
                model.file = textConnection(jags_model), 
                n.chains = 4, 
                n.iter = 20000,
                n.burnin = 10000,
                n.thin = 1,
                DIC = TRUE)

posterior <- samples$BUGSoutput$sims.list$phi
prior <- rbeta(dim(posterior)[1], shape1 = 1, shape2 = 10)

upper_bound <- 0.01

bf <- matrix(NA, nrow = n_participants, ncol = n_groups)

for(g in 1:n_groups){
  for(p in 1:n_participants){
    bf[p,g] <- sum(posterior[,p,g]<=upper_bound)/sum(prior<=upper_bound)
  }
}
colnames(bf) <- c("young", "elderly")
save(bf,file='data/bf-contaminant-model.RData')
