rm(list=ls())
#### Load data and functions ####
load('data/memory.RData')
source("src/Functions.R")
load('data/contaminant_sdt.RData')
original <- sdt
load(file = "data/simulation-multiple-designs-2.Rdata")

#### Participant selection ####
cr <- matrix(NA,nrow=20,ncol=2)
for(aa in 1:2){
  for(ss in 1:20){
    cr[ss,aa] <- (sum(results$res[,ss,aa]==0&results$st[,ss,aa]==1)+
                    sum(results$res[,ss,aa]==1&results$st[,ss,aa]!=1))/192
  }
}

# chose the 3 participants that are closer to the 25, 50 and 75% accuracy in the task from each population
p.sim <- cbind(c(which.min((cr[,1]-quantile(cr[,1],probs = 0.25))^2),
                 which.min((cr[,1]-quantile(cr[,1],probs = 0.5))^2),
                 which.min((cr[,1]-quantile(cr[,1],probs = 0.75))^2)),
               c(which.min((cr[,2]-quantile(cr[,2],probs = 0.25))^2),
                 which.min((cr[,2]-quantile(cr[,2],probs = 0.5))^2),
                 which.min((cr[,2]-quantile(cr[,2],probs = 0.75))^2)))

#### Design choice by participant ####
set.seed(3798731)
# randomly select 15 experimental designs for each participant (total 90)  
designs.p <- array(NA,dim=c(2,3,15,2))
for(aa in 1:2){
  for(ss in 1:3){
    designs.p[aa,ss,,1] <- sample(x = seq(1,20), size = 15, replace = TRUE)
    designs.p[aa,ss,,2] <- rbinom(15,1,0.5)+1
  }
}

#### Distributions for data generation ####
# Posterior mean on after the 192 trials for each participant (used to simulate behavior)
mean.vector <- array(NA,dim=c(7,3,2))
# Variance of the posterior distribution of the parameters scaled down by 1/4
var.matrix <- array(NA,dim=c(7,7,3,2))
for(aa in 1:2){
  for(pp in 1:3){
    mean.vector[,pp,aa] <- original$data$mean[,192,p.sim[pp,aa],aa]    
    var.matrix[,,pp,aa] <- 1/4*original$data$vcov[,,192,p.sim[pp,aa],aa]
  }
}

#### Simulation of behavior for experimental designs ####
sim.data <- list()
sim.data$res <- array(NA,dim=c(192,15,3,2))
sim.data$st <- array(NA,dim=c(192,15,3,2))
for(a in 1:2){
  for(s in 1:3){
    for(d in 1:15){
      for(t in 1:192){
        x <- results$st[t,designs.p[a,s,d,1],designs.p[a,s,d,2]]
        y <- behavior.sim(mean = mean.vector[,s,a],variance = var.matrix[,,s,a],design.point = x)
        sim.data$res[t,d,s,a] <- y
        sim.data$st[t,d,s,a] <- x
      }
    }
  }
}

#### Trial by trial model for simulated observations ####
parameters <- 7
times <- 192
n.sub <- 3
ages <- 2
designs <- 15
sdt$data$mean <- array(NA,dim=c(parameters,times,designs,n.sub,ages))
sdt$data$vcov <- array(NA, dim=c(parameters,parameters,times,designs,n.sub,ages))
sdt$data$ci <- array(NA, dim=c(parameters,6,times,designs,n.sub,ages))
sdt$data$ut <- array(NA,dim=c(times,designs,n.sub,ages))

for(a in 1:ages){
  for(p in 1:3){
    for(d in 1:designs){
      for(t in 1:times){
        y.1 <- sim.data$res[c(seq(1:t)),d,p,a]
        x.1 <- sim.data$st[c(seq(1:t)),d,p,a]
        summaries <- posterior.time.c(rs = y.1,st = x.1)
        sdt$data$mean[,t,d,p,a] <- summaries$posterior$mu
        sdt$data$vcov[,,t,d,p,a] <- summaries$posterior$sigma
        sdt$data$ci[,,t,d,p,a] <- summaries$posterior$ci
      }
      for(t in 1:times){
        sdt$data$ut[t,d,p,a] <- klmnorm(post.mu = as.vector(sdt$data$mean[,t,d,p,a]),
                                        post.sigma = sdt$data$vcov[,,t,d,p,a],
                                        prior.mu = mean.vector[,p,a],
                                        prior.sigma = var.matrix[,,p,a])
      }  
    }
  }
}
