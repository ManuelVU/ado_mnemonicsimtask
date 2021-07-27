# Adaptive design optimization using Contaminant model
# (code is not optimized and may take several days to run!)
# Results are saved in the file 'data/contaminant_sdt.Rdata'
rm(list=ls())
library(R.matlab)

# load functions
source("src/Functions.R")

# load data
load("data/memory.RData")

# number of parameters in the model
parameters <- 7

# number of trials 
times <- 192

# number of participants for each age group
n.sub <- 20

# age groups
ages <- 2

# list to save results from the analisis, data uses the original stimulus and response sequence.
sdt <- list()

# posterior mean by trail for d' and k
sdt$data$mean <- array(NA,dim=c(parameters,times,n.sub,ages))

# posterior covariance matrix by trial for d' and k 
sdt$data$vcov <- array(NA, dim=c(parameters,parameters,times,n.sub,ages))

# posterior credible interval by trial for d' and k
sdt$data$ci <- array(NA, dim=c(parameters,6,times,n.sub,ages))

# posterior mean for of the probability of 'yes' by trial
sdt$pred$mean <- array(NA,dim=c(times,n.sub,ages))

# posterior credible interval of the probability of 'yes' by trial
sdt$pred$ci <- array(NA, dim=c(6,times,n.sub,ages))

# KL divergence between joint posterior distribution of d' and k at trial t and the posterior distribution using all trials
sdt$data$ut <- array(NA,dim=c(times,n.sub,ages))

# posterior probability that trials 1:t are contaminant
sdt$group <- array(NA,dim=c(times,times,n.sub,ages))

for(aa in 1:ages){
  for(pp in 1:n.sub){
    for(tt in 1:times){
      y.1 <- results$res[c(seq(1:tt)),pp,aa]
      x.1 <- results$st[c(seq(1:tt)),pp,aa]
      summaries <- posterior.time.c(rs = y.1,st = x.1,postpred = TRUE)
      sdt$data$mean[,tt,pp,aa] <- summaries$posterior$mu
      sdt$data$vcov[,,tt,pp,aa] <- summaries$posterior$sigma
      sdt$data$ci[,,tt,pp,aa] <- summaries$posterior$ci
      sdt$pred$mean[tt,pp,aa] <- summaries$posterior$pred$mu
      sdt$pred$ci[,tt,pp,aa] <- summaries$posterior$pred$ci
      sdt$group[1:tt,tt,pp,aa] <- summaries$posterior$group
    }
    for(tt in 1:times){
      sdt$data$ut[tt,pp,aa] <- klmnorm(post.mu = as.vector(sdt$data$mean[,tt,pp,aa]),
                                       post.sigma = sdt$data$vcov[,,tt,pp,aa],
                                       prior.mu = sdt$data$mean[,times,pp,aa],
                                       prior.sigma = sdt$data$vcov[,,times,pp,aa])
    }
  }
}

# expected KL between joint posterior of d' and k and the joint prior distribution 
sdt$optimal$kl <- array(NA,dim=c(parameters,times,n.sub,ages))

# these results are the same as in the sdt$data section
sdt$optimal$mean <- array(NA,dim = c(parameters,times,n.sub,ages))
sdt$optimal$vcov <- array(NA,dim=c(parameters,parameters,times,n.sub,ages))
sdt$optimal$ci <- array(NA,dim=c(parameters,6,times,n.sub,ages))
sdt$optimal$ut <- array(NA,dim=c(times,n.sub,ages))

# order in which the original trials where presented in the ADO setting
sdt$optimal$expdes <- array(NA,dim=c(times,n.sub,ages))

# order of the KL divergence values for each stimulus type on each trial
sdt$optimal$storder <- array(NA,dim=c(times,parameters,n.sub,ages))

# posterior probability that trials 1:t are contaminant
sdt$optimal$group <- array(NA,dim=c(times,times,n.sub,ages))
for(aa in 1:ages){
  for(pp in 1:n.sub){
    available.st <- table(results$st[,pp,aa])
    trial.st <- matrix(NA,ncol=parameters,nrow=max(available.st))
    count.tt <- rep(1,7)
    for(kk in 1:parameters){
      trial.st[1:available.st[kk],kk] <- which(results$st[,pp,aa]==kk)
    }
    for(tt in 1:times){
      if(tt==1){
        new.design <- 1
        available.st[results$st[tt,pp,aa]] <- available.st[results$st[tt,pp,aa]]-1
        count.tt[results$st[tt,pp,aa]] <- count.tt[results$st[tt,pp,aa]]+1
        y <- results$res[tt,pp,aa]
        x <- results$st[tt,pp,aa]
        summaries <- kl_sdt.c(rs = y, st = x)
        sdt$optimal$kl[,tt,pp,aa] <- summaries$div
        sdt$optimal$mean[,tt,pp,aa] <- summaries$posterior$mu
        sdt$optimal$vcov[,,tt,pp,aa] <- summaries$posterior$sigma
        sdt$optimal$ci[,,tt,pp,aa] <- summaries$posterior$ci
        sdt$optimal$group[1:tt,tt,pp,aa] <- summaries$posterior$group
      }
      else{
        next.st <- order(sdt$optimal$kl[,(tt-1),pp,aa],decreasing = T)
        bandera <- 0
        q <- 1
        while(bandera==0){
          if(available.st[next.st[q]]==0){
            q <- q+1
            bandera = 0
          }
          else{
            bandera=1
            new.design <- append(new.design,trial.st[count.tt[next.st[q]],next.st[q]])
            available.st[next.st[q]] <- available.st[next.st[q]]-1
            count.tt[next.st[q]] <- count.tt[next.st[q]]+1
          }
        }
      }
      y <- results$res[new.design,pp,aa]
      x <- results$st[new.design,pp,aa]
      summaries <- kl_sdt.c(rs = y, st = x)
      sdt$optimal$kl[,tt,pp,aa] <- summaries$div
      sdt$optimal$mean[,tt,pp,aa] <- summaries$posterior$mu
      sdt$optimal$vcov[,,tt,pp,aa] <- summaries$posterior$sigma
      sdt$optimal$ci[,,tt,pp,aa] <- summaries$posterior$ci
      sdt$optimal$group[1:tt,tt,pp,aa] <- summaries$posterior$group
    }
    sdt$optimal$expdes[,pp,aa] <- new.design
    for(tt in 1:times){
      sdt$optimal$ut[tt,pp,aa] <- klmnorm(post.mu = sdt$optimal$mean[,tt,pp,aa],
                                          post.sigma = sdt$optimal$vcov[,,tt,pp,aa],
                                          prior.mu = sdt$optimal$mean[,times,pp,aa],
                                          prior.sigma = sdt$optimal$vcov[,,times,pp,aa])
      sdt$optimal$storder[tt,,pp,aa] <- order(sdt$optimal$kl[,tt,pp,aa],decreasing = T)
    }
  }
}

# save results of the analysis
save(sdt,file='data/contaminant_sdt.RData')