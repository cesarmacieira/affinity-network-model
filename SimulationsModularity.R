####============================================
#### Aplicação da modularidade - César Macieira
####============================================
rm(list=ls(all=T))

# Packages
if(!require(igraph)){ install.packages("igraph"); require(igraph)}

#Functions
AlocLikelihoodIndC2 = function(n,theta.real,
                               mu1.real,mu2.real,
                               mu1.estimate,mu2.estimate){
  set.seed(13)
  MatrixU = t(cbind(replicate(theta.real[1]*n,sapply(mu1.real,function(x)rbinom(1,1,x))), 
                    replicate(theta.real[2]*n,sapply(mu2.real,function(x)rbinom(1,1,x)))))
  mu.estimate = rbind(mu1.estimate,mu2.estimate)
  MatrixLikelihoods = matrix(NA, ncol = n, nrow = length(theta.real))
  for(k in 1:length(theta.real)){
    for(i in 1:n){
      MatrixLikelihoods[k,i] = prod(mu.estimate[k,]^MatrixU[i,])
    }
  }
  Community1 = Community2 = 0
  Allocation = matrix(nrow = 1, ncol = n) 
  for(i in 1:n){
    if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i]){
      choicedraw = sample(x = c(1,2), size = 1, prob = c(theta.real[1],theta.real[2]))
      if(choicedraw == 1){
        Community1 = Community1 + 1
        Allocation[1,i] = 1
      }else if(choicedraw == 2){
        Community2 = Community2 + 1
        Allocation[1,i] = 2
      }
    }else 
      if(which.max(MatrixLikelihoods[,i]) == 1){
        Community1 = Community1 + 1
        Allocation[1,i] = 1
      }else if(which.max(MatrixLikelihoods[,i]) == 2){
        Community2 = Community2 + 1
        Allocation[1,i] = 2
      }
  }
  results = list()
  results$MatrixLikelihoods = MatrixLikelihoods
  results$Allocation = Allocation
  results$PropAllocation = rbind(Community1/n,Community2/n)
  return(results)
}
AlocLikelihoodIndC3 = function(n,theta.real,
                               mu1.real,mu2.real,mu3.real,
                               mu1.estimate,mu2.estimate,mu3.estimate){
  set.seed(13)
  MatrixU = t(cbind(replicate(theta.real[1]*n,sapply(mu1.real,function(x)rbinom(1,1,x))), 
                    replicate(theta.real[2]*n,sapply(mu2.real,function(x)rbinom(1,1,x))),
                    replicate(theta.real[3]*n,sapply(mu3.real,function(x)rbinom(1,1,x)))))
  mu.estimate = rbind(mu1.estimate,mu2.estimate,mu3.estimate)
  MatrixLikelihoods = matrix(NA, ncol = n, nrow = length(theta.real))
  for(k in 1:length(theta.real)){
    for(i in 1:n){
      MatrixLikelihoods[k,i] = prod(mu.estimate[k,]^MatrixU[i,])
    }
  }
  Community1 = Community2 = Community3 = 0
  Allocation = matrix(nrow = 1, ncol = n) 
  for(i in 1:n){
    if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i]){
      choicedraw = sample(x = c(1,2,3), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[3]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i]){
      choicedraw = sample(x = c(1,2), size = 1, prob = c(theta.real[1],theta.real[2]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1  
      }
      else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i]){
      choicedraw = sample(x = c(1,3), size = 1, prob = c(theta.real[1],theta.real[3]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }
      else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 & MatrixLikelihoods[2,i] == MatrixLikelihoods[3,i]){
      choicedraw = sample(x = c(2,3), size = 1, prob = c(theta.real[2],theta.real[3]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }
      else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1){
      Community1 = Community1 + 1; Allocation[1,i] = 1
    }else if(which.max(MatrixLikelihoods[,i]) == 2){
      Community2 = Community2 + 1; Allocation[1,i] = 2
    }else if(which.max(MatrixLikelihoods[,i]) == 3){
      Community3 = Community3 + 1; Allocation[1,i] = 3
    }
  }
  results = list()
  results$MatrixLikelihoods = MatrixLikelihoods
  results$Allocation = Allocation
  results$PropAllocation = rbind(Community1/n,Community2/n,Community3/n)
  return(results)
  
}
AlocLikelihoodIndC4 = function(n,theta.real,
                               mu1.real,mu2.real,mu3.real,mu4.real,
                               mu1.estimate,mu2.estimate,mu3.estimate,mu4.estimate){
  set.seed(13)
  MatrixU = t(cbind(replicate(theta.real[1]*n,sapply(mu1.real,function(x)rbinom(1,1,x))), 
                    replicate(theta.real[2]*n,sapply(mu2.real,function(x)rbinom(1,1,x))),
                    replicate(theta.real[3]*n,sapply(mu3.real,function(x)rbinom(1,1,x))),
                    replicate(theta.real[4]*n,sapply(mu4.real,function(x)rbinom(1,1,x)))))
  mu.estimate = rbind(mu1.estimate,mu2.estimate,mu3.estimate,mu4.estimate)
  MatrixLikelihoods = matrix(NA, ncol = n, nrow = length(theta.real))
  for(k in 1:length(theta.real)){
    for(i in 1:n){
      MatrixLikelihoods[k,i] = prod(mu.estimate[k,]^MatrixU[i,])
    }
  }
  Community1 = Community2 = Community3 = Community4 = 0
  Allocation = matrix(nrow = 1, ncol = n) 
  for(i in 1:n){
    if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & 
       MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i] &
       MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i]){#Empate entre 1, 2, 3 e 4
      choicedraw = sample(x = c(1,2,3,4), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[3],theta.real[4]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & 
             MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i]){#Empate entre 1, 2 e 3
      choicedraw = sample(x = c(1,2,3), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[3]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i]){#Empate entre 1, 2 e 4
      choicedraw = sample(x = c(1,2,4), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[4]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i]){#Empate entre 1, 3 e 4
      choicedraw = sample(x = c(1,3,4), size = 1, prob = c(theta.real[1],theta.real[3],theta.real[4]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 &
             MatrixLikelihoods[2,i] == MatrixLikelihoods[3,i] &
             MatrixLikelihoods[2,i] == MatrixLikelihoods[4,i]){#Empate entre 2, 3 e 4
      choicedraw = sample(x = c(2,3,4), size = 1, prob = c(theta.real[2],theta.real[3],theta.real[4]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i]){#Empate entre 1 e 2
      choicedraw = sample(x = c(1,2), size = 1, prob = c(theta.real[1],theta.real[2]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1  
      }
      else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i]){#Empate entre 1 e 3
      choicedraw = sample(x = c(1,3), size = 1, prob = c(theta.real[1],theta.real[3]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }
      else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i]){#Empate entre 1 e 4
      choicedraw = sample(x = c(1,4), size = 1, prob = c(theta.real[1],theta.real[4]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }
      else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 & MatrixLikelihoods[2,i] == MatrixLikelihoods[3,i]){#Empate entre 2 e 3
      choicedraw = sample(x = c(2,3), size = 1, prob = c(theta.real[2],theta.real[3]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 & MatrixLikelihoods[2,i] == MatrixLikelihoods[4,i]){#Empate entre 2 e 4
      choicedraw = sample(x = c(2,4), size = 1, prob = c(theta.real[2],theta.real[4]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 3 & MatrixLikelihoods[3,i] == MatrixLikelihoods[4,i]){#Empate entre 3 e 4
      choicedraw = sample(x = c(3,4), size = 1, prob = c(theta.real[3],theta.real[4]))
      if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1){
      Community1 = Community1 + 1; Allocation[1,i] = 1
    }else if(which.max(MatrixLikelihoods[,i]) == 2){
      Community2 = Community2 + 1; Allocation[1,i] = 2
    }else if(which.max(MatrixLikelihoods[,i]) == 3){
      Community3 = Community3 + 1; Allocation[1,i] = 3
    }else if(which.max(MatrixLikelihoods[,i]) == 4){
      Community4 = Community4 + 1; Allocation[1,i] = 4
    }
  }
  results = list()
  results$MatrixLikelihoods = MatrixLikelihoods
  results$Allocation = Allocation
  results$PropAllocation = rbind(Community1/n,Community2/n,Community3/n,Community4/n)
  return(results)
  
}
AlocLikelihoodIndC5 = function(n,theta.real,
                               mu1.real,mu2.real,mu3.real,mu4.real,mu5.real,
                               mu1.estimate,mu2.estimate,mu3.estimate,mu4.estimate,mu5.estimate){
  set.seed(13)
  MatrixU = t(cbind(replicate(theta.real[1]*n,sapply(mu1.real,function(x)rbinom(1,1,x))), 
                    replicate(theta.real[2]*n,sapply(mu2.real,function(x)rbinom(1,1,x))),
                    replicate(theta.real[3]*n,sapply(mu3.real,function(x)rbinom(1,1,x))),
                    replicate(theta.real[4]*n,sapply(mu4.real,function(x)rbinom(1,1,x))),
                    replicate(theta.real[5]*n,sapply(mu5.real,function(x)rbinom(1,1,x)))))
  mu.estimate = rbind(mu1.estimate,mu2.estimate,mu3.estimate,mu4.estimate,mu5.estimate)
  MatrixLikelihoods = matrix(NA, ncol = n, nrow = length(theta.real))
  for(k in 1:length(theta.real)){
    for(i in 1:n){
      MatrixLikelihoods[k,i] = prod(mu.estimate[k,]^MatrixU[i,])
    }
  }
  Community1 = Community2 = Community3 = Community4 = Community5 = 0
  Allocation = matrix(nrow = 1, ncol = n) 
  for(i in 1:n){
    if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & 
       MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i] &
       MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i] &
       MatrixLikelihoods[1,i] == MatrixLikelihoods[5,i]){#Empate entre 1, 2, 3, 4 e 5
      choicedraw = sample(x = c(1,2,3,4,5), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[3],theta.real[4],theta.real[5]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & 
             MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i]){#Empate entre 1, 2, 3 e 4
      choicedraw = sample(x = c(1,2,3,4), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[3],theta.real[4]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & 
             MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[5,i]){#Empate entre 1, 2, 3 e 5
      choicedraw = sample(x = c(1,2,3,5), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[3],theta.real[5]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & 
             MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[5,i]){#Empate entre 1, 2, 4 e 5
      choicedraw = sample(x = c(1,2,4,5), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[4],theta.real[5]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i] & 
             MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[5,i]){#Empate entre 1, 3, 4 e 5
      choicedraw = sample(x = c(1,3,4,5), size = 1, prob = c(theta.real[1],theta.real[3],theta.real[4],theta.real[5]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & 
             MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i]){#Empate entre 1, 2 e 3
      choicedraw = sample(x = c(1,2,3), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[3]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i]){#Empate entre 1, 2 e 4
      choicedraw = sample(x = c(1,2,4), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[4]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[5,i]){#Empate entre 1, 2 e 5
      choicedraw = sample(x = c(1,2,5), size = 1, prob = c(theta.real[1],theta.real[2],theta.real[5]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i]){#Empate entre 1, 3 e 4
      choicedraw = sample(x = c(1,3,4), size = 1, prob = c(theta.real[1],theta.real[3],theta.real[4]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i] &
             MatrixLikelihoods[1,i] == MatrixLikelihoods[5,i]){#Empate entre 1, 3 e 5
      choicedraw = sample(x = c(1,3,5), size = 1, prob = c(theta.real[1],theta.real[3],theta.real[5]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 &
             MatrixLikelihoods[2,i] == MatrixLikelihoods[3,i] &
             MatrixLikelihoods[2,i] == MatrixLikelihoods[4,i]){#Empate entre 2, 3 e 4
      choicedraw = sample(x = c(2,3,4), size = 1, prob = c(theta.real[2],theta.real[3],theta.real[4]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 &
             MatrixLikelihoods[2,i] == MatrixLikelihoods[3,i] &
             MatrixLikelihoods[2,i] == MatrixLikelihoods[5,i]){#Empate entre 2, 3 e 5
      choicedraw = sample(x = c(2,3,5), size = 1, prob = c(theta.real[2],theta.real[3],theta.real[5]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 3 &
             MatrixLikelihoods[3,i] == MatrixLikelihoods[4,i] &
             MatrixLikelihoods[3,i] == MatrixLikelihoods[5,i]){#Empate entre 3, 4 e 5
      choicedraw = sample(x = c(3,4,5), size = 1, prob = c(theta.real[3],theta.real[4],theta.real[5]))
      if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i]){#Empate entre 1 e 2
      choicedraw = sample(x = c(1,2), size = 1, prob = c(theta.real[1],theta.real[2]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1  
      }
      else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i]){#Empate entre 1 e 3
      choicedraw = sample(x = c(1,3), size = 1, prob = c(theta.real[1],theta.real[3]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }
      else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[4,i]){#Empate entre 1 e 4
      choicedraw = sample(x = c(1,4), size = 1, prob = c(theta.real[1],theta.real[4]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }
      else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[5,i]){#Empate entre 1 e 5
      choicedraw = sample(x = c(1,5), size = 1, prob = c(theta.real[1],theta.real[5]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }
      else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 & MatrixLikelihoods[2,i] == MatrixLikelihoods[3,i]){#Empate entre 2 e 3
      choicedraw = sample(x = c(2,3), size = 1, prob = c(theta.real[2],theta.real[3]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 & MatrixLikelihoods[2,i] == MatrixLikelihoods[4,i]){#Empate entre 2 e 4
      choicedraw = sample(x = c(2,4), size = 1, prob = c(theta.real[2],theta.real[4]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 & MatrixLikelihoods[2,i] == MatrixLikelihoods[5,i]){#Empate entre 2 e 5
      choicedraw = sample(x = c(2,5), size = 1, prob = c(theta.real[2],theta.real[5]))
      if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 3 & MatrixLikelihoods[3,i] == MatrixLikelihoods[4,i]){#Empate entre 3 e 4
      choicedraw = sample(x = c(3,4), size = 1, prob = c(theta.real[3],theta.real[4]))
      if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 3 & MatrixLikelihoods[3,i] == MatrixLikelihoods[5,i]){#Empate entre 3 e 5
      choicedraw = sample(x = c(3,5), size = 1, prob = c(theta.real[3],theta.real[5]))
      if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 4 & MatrixLikelihoods[4,i] == MatrixLikelihoods[5,i]){#Empate entre 4 e 5
      choicedraw = sample(x = c(4,5), size = 1, prob = c(theta.real[4],theta.real[5]))
      if(choicedraw == 4){ 
        Community4 = Community4 + 1; Allocation[1,i] = 4 
      }else if(choicedraw == 5){ 
        Community5 = Community5 + 1; Allocation[1,i] = 5 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1){
      Community1 = Community1 + 1; Allocation[1,i] = 1
    }else if(which.max(MatrixLikelihoods[,i]) == 2){
      Community2 = Community2 + 1; Allocation[1,i] = 2
    }else if(which.max(MatrixLikelihoods[,i]) == 3){
      Community3 = Community3 + 1; Allocation[1,i] = 3
    }else if(which.max(MatrixLikelihoods[,i]) == 4){
      Community4 = Community4 + 1; Allocation[1,i] = 4
    }else if(which.max(MatrixLikelihoods[,i]) == 5){
      Community5 = Community5 + 1; Allocation[1,i] = 5
    }
  }
  results = list()
  results$MatrixLikelihoods = MatrixLikelihoods
  results$Allocation = Allocation
  results$PropAllocation = rbind(Community1/n,Community2/n,Community3/n,Community4/n,Community5/n)
  return(results)
  
}

####=======================
#### Simulação 1: Extremos
####=======================
# 1. Inicialização
set.seed(13)
n_pessoas = 10
n_palavras = 9
matriz_medidas = matrix(0, nrow = n_pessoas, ncol = n_palavras)
cluster_real <- rep(NA, n_pessoas)

# 2. Gerar as medidas por grupo
for (i in 1:n_pessoas) {
  if (i <= 3) {
    matriz_medidas[i, 1:3] = runif(3)
    cluster_real[i] = 1
  } else if (i <= 6) {
    matriz_medidas[i, 4:6] = runif(3)
    cluster_real[i] = 2
  } else {
    matriz_medidas[i, 7:9] = runif(3)
    cluster_real[i] = 3
  }
}
matriz_medidas

# 3. Criação da matriz U
matriz_binaria = matrix(0, nrow = n_pessoas, ncol = n_palavras)
for (i in 1:n_pessoas) {
  top_3 = order(matriz_medidas[i, ], decreasing = TRUE)[1:3]
  matriz_binaria[i, top_3] = 1
}
matriz_binaria

# 4. Criar grafo bipartido
linhas = which(matriz_binaria == 1, arr.ind = TRUE)
edges = cbind(paste0("p", linhas[,1]), paste0("w", linhas[,2]))
g_bipartido = graph_from_edgelist(edges, directed = FALSE)
V(g_bipartido)$type = grepl("^w", V(g_bipartido)$name)

# 5. Projetar grafo entre pessoas
g_pessoas = bipartite_projection(g_bipartido)$proj1

# 6. Comunidades e modularidade
comunidades = cluster_louvain(g_pessoas)
modularidade = modularity(comunidades)
cat("Modularidade:", modularidade, "\n")

# 7. Plot do grafo
plot(comunidades, g_pessoas, vertex.label = V(g_pessoas)$name, vertex.color = comunidades$membership,
     vertex.size = 30, edge.width = 2, main = "Grafo de Pessoas com Comunidades (Louvain)")

# 8. Alocações
data.frame(
  id = V(g_pessoas)$name,
  cluster_real = cluster_real[as.numeric(sub("p", "", V(g_pessoas)$name))],
  cluster_modularidade = comunidades$membership
)

AlocTableSim1N1000 = AlocLikelihoodIndC2(n = 60, theta.real = c(0.30,0.70),
                                         mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                         mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                         mu1.estimate = c(Sim1N1000$mu.medium[1,]),
                                         mu2.estimate = c(Sim1N1000$mu.medium[2,]));
AlocTableSim1N1000$PropAllocation

####===================================================================================
#### Simulação 2: uma comunidade com medida alta e duas comunidades com medidas baixas
####===================================================================================
# 1. Inicialização
set.seed(13)
n_pessoas = 100
n_palavras = 90
matriz_medidas = matrix(0, nrow = n_pessoas, ncol = n_palavras)
cluster_real <- rep(NA, n_pessoas)

# 2. Gerar medidas por grupo
for (i in 1:n_pessoas) {
  if (i <= 33) {
    medida1 = runif(n_palavras, min = 0.7, max = 1)  # Grupo 1
    matriz_medidas[i, ] = medida1
    cluster_real[i] = 1
  } else if (i <= 66) {
    medida2 = runif(n_palavras, min = 0, max = 0.2)  # Grupo 2
    matriz_medidas[i, ] = medida2
    cluster_real[i] = 2
  } else {
    medida3 = runif(n_palavras, min = 0.1, max = 0.3)  # Grupo 3
    matriz_medidas[i, ] = medida3
    cluster_real[i] = 3
  }
}
matriz_medidas

# 3. Criar matriz binária U com top 5 palavras por pessoa
matriz_binaria = matrix(0, nrow = n_pessoas, ncol = n_palavras)
for (i in 1:n_pessoas) {
  top_5 = order(matriz_medidas[i, ], decreasing = TRUE)[1:5]
  matriz_binaria[i, top_5] = 1
}
matriz_binaria

# 4. Criar grafo bipartido
linhas = which(matriz_binaria == 1, arr.ind = TRUE)
edges = cbind(paste0("p", linhas[,1]), paste0("w", linhas[,2]))
g_bipartido = graph_from_edgelist(edges, directed = FALSE)
V(g_bipartido)$type = grepl("^w", V(g_bipartido)$name)

# 5. Projetar grafo entre pessoas
g_pessoas = bipartite_projection(g_bipartido)$proj1

# 6. Comunidades e modularidade
comunidades = cluster_louvain(g_pessoas)
modularidade = modularity(comunidades)
cat("Modularidade:", modularidade, "\n")

# 7. Plot do grafo
plot(comunidades, g_pessoas, vertex.label = V(g_pessoas)$name, vertex.color = comunidades$membership,
     vertex.size = 30, edge.width = 2, main = "Grafo de Pessoas com Comunidades (Louvain)")

data.frame(
  id = V(g_pessoas)$name,
  cluster_real = cluster_real[as.numeric(sub("p", "", V(g_pessoas)$name))],
  cluster_modularidade = comunidades$membership
)


set.seed(13)
MatrixU = matriz_binaria
mu.estimate = rbind(medida1,medida2,medida3)
n = n_pessoas
theta.real = c(0.3,0.3,0.4)
MatrixLikelihoods = matrix(NA, ncol = n, nrow = length(theta.real))
for(k in 1:length(theta.real)){
  for(i in 1:n){
    MatrixLikelihoods[k,i] = prod(mu.estimate[k,]^MatrixU[i,])
  }
}
Community1 = Community2 = 0
Allocation = matrix(nrow = 1, ncol = n) 
for(i in 1:n){
  if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i]){
    choicedraw = sample(x = c(1,2), size = 1, prob = c(theta.real[1],theta.real[2]))
    if(choicedraw == 1){
      Community1 = Community1 + 1
      Allocation[1,i] = 1
    }else if(choicedraw == 2){
      Community2 = Community2 + 1
      Allocation[1,i] = 2
    }
  }else 
    if(which.max(MatrixLikelihoods[,i]) == 1){
      Community1 = Community1 + 1
      Allocation[1,i] = 1
    }else if(which.max(MatrixLikelihoods[,i]) == 2){
      Community2 = Community2 + 1
      Allocation[1,i] = 2
    }
}
results = list()
results$MatrixLikelihoods = MatrixLikelihoods
results$Allocation = Allocation
results$PropAllocation = rbind(Community1/n,Community2/n)
