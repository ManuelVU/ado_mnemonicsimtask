# This code simulates the behavior of new participants on the task 
# using the posterior distribution of 3 participants of each age group 
# to generate data, the chosen participants are those whom responses are 
# closest to the 25, 50 and 75% accuracy in the task. (Running this code may take several hours!)
# Results are saved in the file 'data/simulation_multiple_designs.Rdata'

rm(list=ls())
#### Load data and functions ####
load('data/memory.RData')
source("src/od_memory.R")
load('data/contaminant_sdt.RData')
original <- sdt

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
# randomly select 5 experimental designs (no replacement) for each participant (30 total)
designs.p <- array(NA,dim=c(2,3,5,2))
for(aa in 1:2){
  for(ss in 1:3){
    designs.p[aa,ss,,1] <- sample(seq(1,20),5)
    designs.p[aa,ss,,2] <- (rbinom(5,1,0.5)+1)
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
sim.data$res <- array(NA,dim=c(192,5,3,2))
sim.data$st <- array(NA,dim=c(192,5,3,2))
for(a in 1:2){
  for(s in 1:3){
    for(d in 1:5){
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
designs <- 5
sdt <- list()
sdt$data$mean <- array(NA,dim=c(parameters,times,designs,n.sub,ages))
sdt$data$vcov <- array(NA, dim=c(parameters,parameters,times,designs,n.sub,ages))
sdt$data$ci <- array(NA, dim=c(parameters,6,times,designs,n.sub,ages))
sdt$data$ut <- array(NA,dim=c(times,designs,n.sub,ages))

for(a in 1:ages){
  for(p in 1:3){
    for(d in 1:5){
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

#### Trial by trial optimal design ####
parameters <- 7
times <- 192
n.sub <- 3
ages <- 2
sdt$optimal$kl <- array(NA,dim=c(parameters,times,n.sub,ages))
sdt$optimal$mean <- array(NA,dim = c(parameters,times,n.sub,ages))
sdt$optimal$vcov <- array(NA,dim=c(parameters,parameters,times,n.sub,ages))
sdt$optimal$ci <- array(NA,dim=c(parameters,6,times,n.sub,ages))
sdt$optimal$ut <- array(NA,dim=c(times,n.sub,ages))
sdt$optimal$expdes <- array(NA,dim=c(times,n.sub,ages))
sdt$optimal$storder <- array(NA,dim=c(times,parameters,n.sub,ages))

for(aa in 1:ages){
  for(pp in 1:n.sub){
    for(tt in 1:times){
      if(tt==1){
        new.design <- sample(seq(1,7),size = 1)
        dp <- new.design
        obs <- behavior.sim(mean = mean.vector[,pp,aa], variance = var.matrix[,,pp,aa],design.point = dp)
        summaries <- kl_sdt.c(rs = obs, st = dp)
        sdt$optimal$kl[,tt,pp,aa] <- summaries$div
        sdt$optimal$mean[,tt,pp,aa] <- summaries$posterior$mu
        sdt$optimal$vcov[,,tt,pp,aa] <- summaries$posterior$sigma
        sdt$optimal$ci[,,tt,pp,aa] <- summaries$posterior$ci
      }
      else{
        next.st <- order(sdt$optimal$kl[,(tt-1),pp,aa],decreasing = T)
        new.design <- append(new.design,next.st[1])
        dp <- next.st[1]
        y.1  <- behavior.sim(mean = mean.vector[,pp,aa], variance = var.matrix[,,pp,aa], design.point = dp)
        obs <- append(obs,y.1)
        summaries <- kl_sdt.c(rs = obs, st = new.design)
        sdt$optimal$kl[,tt,pp,aa] <- summaries$div
        sdt$optimal$mean[,tt,pp,aa] <- summaries$posterior$mu
        sdt$optimal$vcov[,,tt,pp,aa] <- summaries$posterior$sigma
        sdt$optimal$ci[,,tt,pp,aa] <- summaries$posterior$ci
      }
    }
    sdt$optimal$expdes[,pp,aa] <- new.design
    for(tt in 1:192){
      sdt$optimal$ut[tt,pp,aa] <- klmnorm(post.mu = as.vector(sdt$optimal$mean[,tt,pp,aa]),
                                          post.sigma = sdt$optimal$vcov[,,tt,pp,aa],
                                          prior.mu = mean.vector[,pp,aa],
                                          prior.sigma = var.matrix[,,pp,aa])
    }
  }
}

# save results 
save(sdt, file = 'data/simulation_multiple_designs.Rdata')
