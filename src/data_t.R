# This file reads and re arrenges the data from the mnemonic similarity task.
# The result from this file is saved in 'data/memory.Rdata'. This file is used 
# as input for all functions and analysis that require experimental data.
rm(list=ls())
library(R.matlab)
# Read original data file
star <- readMat('data/StarkExp1.mat')
# Format for memory array: 1: 192 trials, 2: 6 variables 
# (trial type, lure bin, response(1=new), accuracy, response time, confidence)
# 3: 20 participants, 4: 2 age groups as young-old.
memory <- array(NA, dim=c(192,6,20,2))
for(a in 1:2){
  for(i in 1:(dim(star$d[[1]])[1])/2){
    ac <- a-1
    memory[,1,i,a] <- star$d[[3]][(i+20*ac),] # trial type
    memory[,2,i,a] <- star$d[[4]][(i+20*ac),] # lure bin
    memory[,3,i,a] <- star$d[[5]][(i+20*ac),] # response (1 = Old)
    memory[,4,i,a] <- star$d[[6]][(i+20*ac),] # accuracy (1 = Correct rejection)
    memory[,5,i,a] <- star$d[[7]][(i+20*ac),] # response time
    memory[,6,i,a] <- star$d[[8]][(i+20*ac),] # confidence
  }
}
# coding for response (var 3, 1=new 0=old)
# lure bin for new items changed from 0 to 6
for(a in 1:2){
  for(i in 1:20){
    memory[,3,i,a] <- (memory[,3,i,a]-1)
    memory[which(memory[,1,i,a]==2),2,i,a] <- rep(6,64)
  }
}
memory[,2,,] <-memory[,2,,]+1
y <- array(NA,dim=c(192,20,2))
stimulus <- array(NA, dim=c(192,20,2))
for(a in 1:2){
  for(i in 1:20){
    y[,i,a] <- memory[,3,i,a]
    stimulus[,i,a] <- memory[,2,i,a]
  }
}
results <- list()
results$res <- y
results$st <- stimulus
save(results, file='data/memory.RData')