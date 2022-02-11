load(here::here("data/contaminant-mean-theta.RData"))
load(here::here("data/memory.RData"))
age_groups <- dim(theta_hat)[1]
n_participants <- dim(theta_hat)[2]

boundary <- seq(0,1,0.18)
prop_response <- matrix(data = NA, nrow = age_groups, ncol = length(boundary)-1)
std_error <- matrix(data = NA, nrow = age_groups, ncol = length(boundary)-1)

for(a in 1:age_groups){
  for(i in 2:length(boundary)){
    size_tmp <- c()
    new_tmp <- c()
    for(p in 1:n_participants){
      size_tmp <- append(x = size_tmp, 
                         values = length(which(theta_hat[a,p,] >= boundary[i-1] &  
                                               theta_hat[a,p,] < boundary[i])))
      new_tmp <- append(x = new_tmp, 
                        values = sum(results$res[which(theta_hat[a,p,] >= boundary[i-1] &  
                                           theta_hat[a,p,] < boundary[i]),p,a]))
      
      
    }
    prop_response[a,i-1] <- sum(new_tmp)/sum(size_tmp)
    std_error[a,i-1] <- sqrt((prop_response[a,i-1]*(1-prop_response[a,i-1]))/sum(size_tmp))
  }
}

prop_response <- ifelse(test = is.finite(prop_response), yes = prop_response,
                        no = 0)

accuracy <- array(NA, dim = c(192,20,2))

for(a in 1:age_groups){
  for(p in 1:n_participants){
    accuracy[,p,a] <- ifelse(test = theta_hat[a,p,] > 0.5 & results$res[,p,a] == 1, yes = 1,
                             no = ifelse(
                               test = theta_hat[a,p,] <= 0.5 & results$res[,p,a] == 0, 
                               yes = 1, no = 0))
  }
}

accuracy_prop <- rbind(colMeans(accuracy))

save(accuracy_prop,file='data/contaminant-accuracy.RData')