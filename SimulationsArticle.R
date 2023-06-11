####===============================
#### EM algorithm - César Macieira
####===============================
####==============================================================
#### Simulations of several Scenarios for the initial values of µ
#### µ1 = [0.3 0.4 0.1 0.7 0.0] e µ2 = [0.4 0.0 0.0 0.1 0.3]
#### θ1 = 0.3 e θ2 = 0.7
####==============================================================
rm(list=ls(all=T))

#Packages
if(!require(openxlsx)){ install.packages("openxlsx"); require(openxlsx)}
if(!require(ggplot2)){ install.packages("ggplot2"); require(ggplot2)}
if(!require(tidyverse)){ install.packages("tidyverse"); require(tidyverse)}

# Parameters
# m = number of dictionary words
# µ1 = probability measure for community 1 with size m
# µ2 = probability measure for community 2 with size m
# θ1 = proportion of people in the community 1
# θ2 = proportion of people in the community 2
# θ1 + θ2 = 1
# n1 = number of individuals who choose with µ1
# n2 = number of individuals who choose with µ2
# n = n1 + n2

EMalgorithm = function(MatrixU,C,tol = 1e-8,mu.initialvalues){
  m = ncol(MatrixU)
  n = nrow(MatrixU)
  mu.estimate = mu.initialvalues
  theta.estimate = rep(1/C,C)
  iterations = 0
  Estart = 1
  E = 0 
  
  while(abs(Estart - E) >= tol){ 
    mue = mu.estimate
    thetae = theta.estimate
    Estart = E
    mumtx = list()
    for(i in 1:C){
      mumtx[[i]] = t(replicate(n,mu.estimate[i,]))
    }
    p1 = list()
    p2 = matrix(nrow = n, ncol = C)
    for(k in 1:C){ # Calculate the product of the numerator of Tik
      p1[[k]] = (mumtx[[k]]^MatrixU) * ((1-mumtx[[k]])^(1-MatrixU))
      for(i in 1:n){
        p2[i,k] = prod(p1[[k]][i,])
      }
    }
    ptheta = matrix(rep(theta.estimate,n),ncol = C,byrow = T)
    p3 = p2*ptheta # Calculate the numerator of Tik
    #Calculate Tik
    Tik = p3/matrix(rep(rowSums(p3),C),ncol = C)
    # Calculate Q(µ|µ(t))
    part1 = sum(Tik*log(ptheta))
    part2 = part3 = 0 
    for(k in 1:C){
      part2 = part2 + sum((Tik*rowSums(log(mumtx[[k]]+1e-300)*MatrixU))[,k])
      part3 = part3 + sum((Tik*rowSums(log(1-mumtx[[k]]+1e-300)*(1-MatrixU)))[,k])
    }
    E = part1 + part2 + part3
    #Estimate θ
    theta.estimate = colMeans(Tik)
    #Estimate µ
    Tikmtx = list()
    for(k in 1:C){
      Tikmtx[[k]] = replicate(m,Tik[,k])
    }
    mu.estimate = matrix(nrow = C, ncol = m)
    for(k in 1:C){
      mu.estimate[k,] = t(colSums(Tikmtx[[k]]*MatrixU))/colSums(Tikmtx[[k]])
    }
    iterations = iterations + 1
  }
  output = list()
  output$iterations = iterations
  output$logLikelihood = E
  output$theta = thetae
  output$mu = round(mue,3)
  return(output)
}
EMsimulations = function(n,C,nsim,theta.real,mu1.real,mu2.real,mu1.initial,mu2.initial){
  mu.true = rbind(mu1.real,mu2.real)
  theta.sim = matrix(ncol = C, nrow = nsim)
  mu.sim = list()
  Likelihood = matrix(ncol = 1, nrow = nsim)
  for(sim in 1:nsim){ 
    MatrixU = t(cbind(replicate(theta.real[1]*n,sapply(mu1.real,function(x)rbinom(1,1,x))), 
                      replicate(theta.real[2]*n,sapply(mu2.real,function(x)rbinom(1,1,x)))))
    mu.initialvalues = rbind(mu1.initial,mu2.initial)
    EMestimate = EMalgorithm(MatrixU,C,mu.initialvalues=rbind(mu1.initial,mu2.initial))
    Likelihood[sim,1] = EMestimate$logLikelihood
    for(k in 1:C){
      theta.sim[sim,k] = EMestimate$theta[k]
      mu.sim[[sim]] = EMestimate$mu
    }
  }
  mu.medium = matrix(nrow = C, ncol = ncol(MatrixU))
  stdev = matrix(nrow = C, ncol = ncol(MatrixU))
  for(i in 1:C){
    for(j in 1:ncol(MatrixU)){
      mu.sum = 0
      auxiliar = matrix(nrow = 1, ncol = nsim)
      for(sim in 1:nsim){
        mu.sum = mu.sum + mu.sim[[sim]][i,j]
        auxiliar[1,sim] = mu.sim[[sim]][i,j]
      }
      mu.medium[i,j] = mu.sum/nsim
      stdev[i,j] = sd(auxiliar)
    }
  }
  
  parameters.est = list()
  parameters.est$MatrixU = MatrixU
  parameters.est$Likelihood = Likelihood
  parameters.est$theta = theta.sim
  parameters.est$mu = mu.sim
  parameters.est$mu.medium = mu.medium
  parameters.est$mu.stdev = stdev
  return(parameters.est)
  
}
EMsimulations3Communities = function(n,C,nsim,theta.real,mu1.real,mu2.real,mu3.real,
                                     mu1.initial,mu2.initial,mu3.initial){
  mu.true = rbind(mu1.real,mu2.real,mu3.real)
  theta.sim = matrix(ncol = C, nrow = nsim)
  mu.sim = list()
  for(sim in 1:nsim){ 
    MatrixU = t(cbind(replicate(theta.real[1]*n,sapply(mu1.real,function(x)rbinom(1,1,x))), 
                      replicate(theta.real[2]*n,sapply(mu2.real,function(x)rbinom(1,1,x))),
                      replicate(theta.real[3]*n,sapply(mu3.real,function(x)rbinom(1,1,x)))))
    mu.initialvalues = rbind(mu1.initial,mu2.initial,mu3.initial)
    EMestimate = EMalgorithm(MatrixU,C,mu.initialvalues=rbind(mu1.initial,mu2.initial,mu3.initial))
    for(k in 1:C){
      theta.sim[sim,k] = EMestimate$theta[k]
      mu.sim[[sim]] = EMestimate$mu
    }
  }
  mu.medium = matrix(nrow = C, ncol = ncol(MatrixU))
  stdev = matrix(nrow = C, ncol = ncol(MatrixU))
  for(i in 1:C){
    for(j in 1:ncol(MatrixU)){
      mu.sum = 0
      auxiliar = matrix(nrow = 1, ncol = nsim)
      for(sim in 1:nsim){
        mu.sum = mu.sum + mu.sim[[sim]][i,j]
        auxiliar[1,sim] = mu.sim[[sim]][i,j]
      }
      mu.medium[i,j] = mu.sum/nsim
      stdev[i,j] = sd(auxiliar)
    }
  }
  
  parameters.est = list()
  parameters.est$MatrixU = MatrixU
  parameters.est$theta = theta.sim
  parameters.est$mu = mu.sim
  parameters.est$mu.medium = mu.medium
  parameters.est$mu.stdev = stdev
  return(parameters.est)
}
EMsimulations4Communities = function(n,C,nsim,theta.real,mu1.real,mu2.real,mu3.real,mu4.real,
                                     mu1.initial,mu2.initial,mu3.initial,mu4.initial){
  mu.true = rbind(mu1.real,mu2.real,mu3.real,mu4.real)
  theta.sim = matrix(ncol = C, nrow = nsim)
  mu.sim = list()
  for(sim in 1:nsim){ 
    MatrixU = t(cbind(replicate(theta.real[1]*n,sapply(mu1.real,function(x)rbinom(1,1,x))), 
                      replicate(theta.real[2]*n,sapply(mu2.real,function(x)rbinom(1,1,x))),
                      replicate(theta.real[3]*n,sapply(mu3.real,function(x)rbinom(1,1,x))),
                      replicate(theta.real[4]*n,sapply(mu4.real,function(x)rbinom(1,1,x)))))
    mu.initialvalues = rbind(mu1.initial,mu2.initial,mu3.initial,mu4.initial)
    EMestimate = EMalgorithm(MatrixU,C,mu.initialvalues=rbind(mu1.initial,mu2.initial,
                                                              mu3.initial,mu4.initial))
    for(k in 1:C){
      theta.sim[sim,k] = EMestimate$theta[k]
      mu.sim[[sim]] = EMestimate$mu
    }
  }
  mu.medium = matrix(nrow = C, ncol = ncol(MatrixU))
  stdev = matrix(nrow = C, ncol = ncol(MatrixU))
  for(i in 1:C){
    for(j in 1:ncol(MatrixU)){
      mu.sum = 0
      auxiliar = matrix(nrow = 1, ncol = nsim)
      for(sim in 1:nsim){
        mu.sum = mu.sum + mu.sim[[sim]][i,j]
        auxiliar[1,sim] = mu.sim[[sim]][i,j]
      }
      mu.medium[i,j] = mu.sum/nsim
      stdev[i,j] = sd(auxiliar)
    }
  }
  
  parameters.est = list()
  parameters.est$MatrixU = MatrixU
  parameters.est$theta = theta.sim
  parameters.est$mu = mu.sim
  parameters.est$mu.medium = mu.medium
  parameters.est$mu.stdev = stdev
  return(parameters.est)
}
EMsimulations5Communities = function(n,C,nsim,theta.real,
                                     mu1.real,mu2.real,mu3.real,mu4.real,mu5.real,
                                     mu1.initial,mu2.initial,mu3.initial,mu4.initial,mu5.initial){
  mu.true = rbind(mu1.real,mu2.real,mu3.real,mu4.real,mu5.real)
  theta.sim = matrix(ncol = C, nrow = nsim)
  mu.sim = list()
  for(sim in 1:nsim){ 
    MatrixU = t(cbind(replicate(theta.real[1]*n,sapply(mu1.real,function(x)rbinom(1,1,x))), 
                      replicate(theta.real[2]*n,sapply(mu2.real,function(x)rbinom(1,1,x))),
                      replicate(theta.real[3]*n,sapply(mu3.real,function(x)rbinom(1,1,x))),
                      replicate(theta.real[4]*n,sapply(mu4.real,function(x)rbinom(1,1,x))),
                      replicate(theta.real[5]*n,sapply(mu5.real,function(x)rbinom(1,1,x)))))
    mu.initialvalues = rbind(mu1.initial,mu2.initial,mu3.initial,mu4.initial,mu5.initial)
    EMestimate = EMalgorithm(MatrixU,C,mu.initialvalues=rbind(mu1.initial,mu2.initial,
                                                              mu3.initial,mu4.initial,
                                                              mu5.initial))
    for(k in 1:C){
      theta.sim[sim,k] = EMestimate$theta[k]
      mu.sim[[sim]] = EMestimate$mu
    }
  }
  mu.medium = matrix(nrow = C, ncol = ncol(MatrixU))
  stdev = matrix(nrow = C, ncol = ncol(MatrixU))
  for(i in 1:C){
    for(j in 1:ncol(MatrixU)){
      mu.sum = 0
      auxiliar = matrix(nrow = 1, ncol = nsim)
      for(sim in 1:nsim){
        mu.sum = mu.sum + mu.sim[[sim]][i,j]
        auxiliar[1,sim] = mu.sim[[sim]][i,j]
      }
      mu.medium[i,j] = mu.sum/nsim
      stdev[i,j] = sd(auxiliar)
    }
  }
  
  parameters.est = list()
  parameters.est$MatrixU = MatrixU
  parameters.est$theta = theta.sim
  parameters.est$mu = mu.sim
  parameters.est$mu.medium = mu.medium
  parameters.est$mu.stdev = stdev
  return(parameters.est)
  
}
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

####=====================================
#### µ - set.seed(13)
#### Change of initial values
#### µ1 = [0.00, 0.10, 0.40, 0.10, 0.00]
#### µ2 = [0.30, 0.40, 0.10, 0.90, 0.80]
####=====================================
####==========================================================
#### Scenario 1 - Number of simulations = 1000 - set.seed(10)
####==========================================================
set.seed(13);
Sim1N1000 = EMsimulations(n = 60, C = 2, nsim = 1000,
                          theta.real = c(0.3,0.7),
                          mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                          mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                          mu1.initial = c(0.1,0.2,0.5,0.1,0.2),
                          mu2.initial = c(0.1,0.3,0.2,0.8,0.6));
AlocTableSim1N1000 = AlocLikelihoodIndC2(n = 60, theta.real = c(0.30,0.70),
                                         mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                         mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                         mu1.estimate = c(Sim1N1000$mu.medium[1,]),
                                         mu2.estimate = c(Sim1N1000$mu.medium[2,]));
TableSim1N1000 = cbind(rbind(mean(Sim1N1000$theta[,1]),
                             mean(Sim1N1000$theta[,2])),
                       Sim1N1000$mu.medium,
                       rbind(sd(Sim1N1000$theta[,1]),sd(Sim1N1000$theta[,2])),
                       Sim1N1000$mu.stdev,AlocTableSim1N1000$PropAllocation);TableSim1N1000

####==========================================================
#### Scenario 2 - Number of simulations = 1000 - set.seed(13)
####==========================================================
set.seed(13)
Sim2N1000 = EMsimulations(n = 60, C = 2, nsim = 20,
                          theta.real = c(0.3,0.7),
                          mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                          mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                          mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                          mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocTableSim2N1000 = AlocLikelihoodIndC2(n = 60, theta.real = c(0.30,0.70),
                                         mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                         mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                         mu1.estimate = c(Sim2N1000$mu.medium[1,]),
                                         mu2.estimate = c(Sim2N1000$mu.medium[2,]));
TableSim2N1000 = cbind(rbind(mean(Sim2N1000$theta[,1]),
                             mean(Sim2N1000$theta[,2])),
                       Sim2N1000$mu.medium,
                       rbind(sd(Sim2N1000$theta[,1]),sd(Sim2N1000$theta[,2])),
                       Sim2N1000$mu.stdev,AlocTableSim2N1000$PropAllocation);TableSim2N1000

####==========================================================
#### Scenario 3 - Number of simulations = 1000 - set.seed(13)
####==========================================================
set.seed(13)
Sim3N1000 = EMsimulations(n = 60, C = 2, nsim = 1000,
                          theta.real = c(0.3,0.7),
                          mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                          mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                          mu1.initial = c(0.2,0.2,0.2,0.2,0.2),
                          mu2.initial = c(0.5,0.5,0.5,0.5,0.5));
AlocTableSim3N1000 = AlocLikelihoodIndC2(n = 60, theta.real = c(0.30,0.70),
                                         mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                         mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                         mu1.estimate = c(Sim3N1000$mu.medium[1,]),
                                         mu2.estimate = c(Sim3N1000$mu.medium[2,]));
TableSim3N1000 = cbind(rbind(mean(Sim3N1000$theta[,1]),
                             mean(Sim3N1000$theta[,2])),
                       Sim3N1000$mu.medium,
                       rbind(sd(Sim3N1000$theta[,1]),sd(Sim3N1000$theta[,2])),
                       Sim3N1000$mu.stdev,AlocTableSim3N1000$PropAllocation);TableSim3N1000

####==========================================================
#### Scenario 4 - Number of simulations = 1000 - set.seed(13)
####==========================================================
set.seed(13)
Sim4N1000 = EMsimulations(n = 60, C = 2, nsim = 1000,
                          theta.real = c(0.3,0.7),
                          mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                          mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                          mu1.initial = c(0.1,0.5,0.2,0.4,0.2),
                          mu2.initial = c(0.2,0.3,0.3,0.5,0.6));
AlocTableSim4N1000 = AlocLikelihoodIndC2(n = 60, theta.real = c(0.30,0.70),
                                         mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                         mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                         mu1.estimate = c(Sim4N1000$mu.medium[1,]),
                                         mu2.estimate = c(Sim4N1000$mu.medium[2,]));
TableSim4N1000 = cbind(rbind(mean(Sim4N1000$theta[,1]),
                             mean(Sim4N1000$theta[,2])),
                       Sim4N1000$mu.medium,
                       rbind(sd(Sim4N1000$theta[,1]),sd(Sim4N1000$theta[,2])),
                       Sim4N1000$mu.stdev,AlocTableSim4N1000$PropAllocation);TableSim4N1000

####==========================================================
#### Scenario 5 - Number of simulations = 1000 - set.seed(13)
####==========================================================
set.seed(13)
Sim5N1000 = EMsimulations(n = 60, C = 2, nsim = 1000,
                          theta.real = c(0.3,0.7),
                          mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                          mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                          mu1.initial = c(0.1,0.3,0.1,0.4,0.7),
                          mu2.initial = c(0.4,0.1,0.3,0.8,0.9));
AlocTableSim5N1000 = AlocLikelihoodIndC2(n = 60, theta.real = c(0.30,0.70),
                                         mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                         mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                         mu1.estimate = c(Sim5N1000$mu.medium[1,]),
                                         mu2.estimate = c(Sim5N1000$mu.medium[2,]));
TableSim5N1000 = cbind(rbind(mean(Sim5N1000$theta[,1]),
                             mean(Sim5N1000$theta[,2])),
                       Sim5N1000$mu.medium,
                       rbind(sd(Sim5N1000$theta[,1]),sd(Sim5N1000$theta[,2])),
                       Sim5N1000$mu.stdev,AlocTableSim5N1000$PropAllocation);TableSim5N1000

####==========================================================
#### Scenario 6 - Number of simulations = 1000 - set.seed(13)
####==========================================================
set.seed(13)
Sim6N1000 = EMsimulations(n = 60, C = 2, nsim = 1000,
                          theta.real = c(0.3,0.7),
                          mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                          mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                          mu1.initial = c(0.1,0.2,0.1,0.2,0.1),
                          mu2.initial = c(0.1,0.7,0.2,0.2,0.1));
AlocTableSim6N1000 = AlocLikelihoodIndC2(n = 60, theta.real = c(0.30,0.70),
                                         mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                         mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                         mu1.estimate = c(Sim6N1000$mu.medium[1,]),
                                         mu2.estimate = c(Sim6N1000$mu.medium[2,]));
TableSim6N1000 = cbind(rbind(mean(Sim6N1000$theta[,1]),
                             mean(Sim6N1000$theta[,2])),
                       Sim6N1000$mu.medium,
                       rbind(sd(Sim6N1000$theta[,1]),sd(Sim6N1000$theta[,2])),
                       Sim6N1000$mu.stdev,AlocTableSim6N1000$PropAllocation);TableSim6N1000

####============================
#### µ - Construction of Tables
####============================
Table7_appendix = rbind(TableSim1N1000,TableSim2N1000,TableSim3N1000,
                        TableSim4N1000,TableSim5N1000,TableSim6N1000)

####=====================================
#### N - set.seed(13)
#### Changing the number of vertices
#### µ1 = [0.00, 0.10, 0.40, 0.10, 0.00]
#### µ2 = [0.30, 0.40, 0.10, 0.90, 0.80]
####=====================================
####=====================================================
#### Scenario 1 - Number of vertices = 20 - set.seed(13)
####=====================================================
# 20 vértices 
set.seed(13)
Scenario1Sim1000N20 = EMsimulations(n = 20, C = 2, nsim = 1000,
                                    theta.real = c(0.3,0.7),
                                    mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                    mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                    mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                                    mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocScenario1Sim1000N20 = AlocLikelihoodIndC2(n = 20, theta.real = c(0.30,0.70),
                                              mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                              mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                              mu1.estimate = c(Scenario1Sim1000N20$mu.medium[1,]),
                                              mu2.estimate = c(Scenario1Sim1000N20$mu.medium[2,]));
TableScenario1Sim1000N20 = cbind(rbind(mean(Scenario1Sim1000N20$theta[,1]),
                                       mean(Scenario1Sim1000N20$theta[,2])),
                                 Scenario1Sim1000N20$mu.medium,
                                 rbind(sd(Scenario1Sim1000N20$theta[,1]),sd(Scenario1Sim1000N20$theta[,2])),
                                 Scenario1Sim1000N20$mu.stdev,AlocScenario1Sim1000N20$PropAllocation);TableScenario1Sim1000N20

####=====================================================
#### Scenario 2 - Number of vertices = 60 - set.seed(13)
####=====================================================
set.seed(13)
Scenario2Sim1000N60 = EMsimulations(n = 60, C = 2, nsim = 1000,
                                    theta.real = c(0.3,0.7),
                                    mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                    mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                    mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                                    mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocScenario2Sim1000N60 = AlocLikelihoodIndC2(n = 60, theta.real = c(0.30,0.70),
                                              mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                              mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                              mu1.estimate = c(Scenario2Sim1000N60$mu.medium[1,]),
                                              mu2.estimate = c(Scenario2Sim1000N60$mu.medium[2,]));
TableScenario2Sim1000N60 = cbind(rbind(mean(Scenario2Sim1000N60$theta[,1]),
                                       mean(Scenario2Sim1000N60$theta[,2])),
                                 Scenario2Sim1000N60$mu.medium,
                                 rbind(sd(Scenario2Sim1000N60$theta[,1]),sd(Scenario2Sim1000N60$theta[,2])),
                                 Scenario2Sim1000N60$mu.stdev,AlocScenario2Sim1000N60$PropAllocation);TableScenario2Sim1000N60

####======================================================
#### Scenario 3 - Number of vertices = 100 - set.seed(13)
####======================================================
set.seed(13)
Scenario3Sim1000N100 = EMsimulations(n = 100, C = 2, nsim = 1000,
                                     theta.real = c(0.3,0.7),
                                     mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                     mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                     mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                                     mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocScenario3Sim1000N100 = AlocLikelihoodIndC2(n = 100, theta.real = c(0.30,0.70),
                                               mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                               mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                               mu1.estimate = c(Scenario3Sim1000N100$mu.medium[1,]),
                                               mu2.estimate = c(Scenario3Sim1000N100$mu.medium[2,]));
TableScenario3Sim1000N100 = cbind(rbind(mean(Scenario3Sim1000N100$theta[,1]),
                                        mean(Scenario3Sim1000N100$theta[,2])),
                                  Scenario3Sim1000N100$mu.medium,
                                  rbind(sd(Scenario3Sim1000N100$theta[,1]),sd(Scenario3Sim1000N100$theta[,2])),
                                  Scenario3Sim1000N100$mu.stdev,AlocScenario3Sim1000N100$PropAllocation);TableScenario3Sim1000N100

####======================================================
#### Scenario 4 - Number of vertices = 150 - set.seed(13)
####======================================================
set.seed(13)
Scenario4Sim1000N150 = EMsimulations(n = 150, C = 2, nsim = 1000,
                                     theta.real = c(0.3,0.7),
                                     mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                     mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                     mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                                     mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocScenario4Sim1000N150 = AlocLikelihoodIndC2(n = 150, theta.real = c(0.30,0.70),
                                               mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                               mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                               mu1.estimate = c(Scenario4Sim1000N150$mu.medium[1,]),
                                               mu2.estimate = c(Scenario4Sim1000N150$mu.medium[2,]));
TableScenario4Sim1000N150 = cbind(rbind(mean(Scenario4Sim1000N150$theta[,1]),
                                        mean(Scenario4Sim1000N150$theta[,2])),
                                  Scenario4Sim1000N150$mu.medium,
                                  rbind(sd(Scenario4Sim1000N150$theta[,1]),sd(Scenario4Sim1000N150$theta[,2])),
                                  Scenario4Sim1000N150$mu.stdev,AlocScenario4Sim1000N150$PropAllocation);TableScenario4Sim1000N150

####======================================================
#### Scenario 5 - Number of vertices = 200 - set.seed(13)
####======================================================
set.seed(13)
Scenario5Sim1000N200 = EMsimulations(n = 200, C = 2, nsim = 1000,
                                     theta.real = c(0.3,0.7),
                                     mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                     mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                     mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                                     mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocScenario5Sim1000N200 = AlocLikelihoodIndC2(n = 200, theta.real = c(0.30,0.70),
                                               mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                               mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                               mu1.estimate = c(Scenario5Sim1000N200$mu.medium[1,]),
                                               mu2.estimate = c(Scenario5Sim1000N200$mu.medium[2,]));
TableScenario5Sim1000N200 = cbind(rbind(mean(Scenario5Sim1000N200$theta[,1]),
                                        mean(Scenario5Sim1000N200$theta[,2])),
                                  Scenario5Sim1000N200$mu.medium,
                                  rbind(sd(Scenario5Sim1000N200$theta[,1]),sd(Scenario5Sim1000N200$theta[,2])),
                                  Scenario5Sim1000N200$mu.stdev,AlocScenario5Sim1000N200$PropAllocation);TableScenario5Sim1000N200

####======================================================
#### Scenario 6 - Number of vertices = 500 - set.seed(13)
####======================================================
set.seed(13)
Scenario6Sim1000N500 = EMsimulations(n = 500, C = 2, nsim = 1000,
                                     theta.real = c(0.3,0.7),
                                     mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                     mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                     mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                                     mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocScenario6Sim1000N500 = AlocLikelihoodIndC2(n = 500, theta.real = c(0.30,0.70),
                                               mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                               mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                               mu1.estimate = c(Scenario6Sim1000N500$mu.medium[1,]),
                                               mu2.estimate = c(Scenario6Sim1000N500$mu.medium[2,]));
TableScenario6Sim1000N500 = cbind(rbind(mean(Scenario6Sim1000N500$theta[,1]),
                                        mean(Scenario6Sim1000N500$theta[,2])),
                                  Scenario6Sim1000N500$mu.medium,
                                  rbind(sd(Scenario6Sim1000N500$theta[,1]),sd(Scenario6Sim1000N500$theta[,2])),
                                  Scenario6Sim1000N500$mu.stdev,AlocScenario6Sim1000N500$PropAllocation);TableScenario6Sim1000N500

####=======================================================
#### Scenario 7 - Number of vertices = 1000 - set.seed(13)
####=======================================================
set.seed(13)
Scenario7Sim1000N1000 = EMsimulations(n = 1000, C = 2, nsim = 1000,
                                      theta.real = c(0.3,0.7),
                                      mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                      mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                      mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                                      mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocScenario7Sim1000N1000 = AlocLikelihoodIndC2(n = 1000, theta.real = c(0.30,0.70),
                                                mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                                mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                                mu1.estimate = c(Scenario7Sim1000N1000$mu.medium[1,]),
                                                mu2.estimate = c(Scenario7Sim1000N1000$mu.medium[2,]));
TableScenario7Sim1000N1000 = cbind(rbind(mean(Scenario7Sim1000N1000$theta[,1]),
                                         mean(Scenario7Sim1000N1000$theta[,2])),
                                   Scenario7Sim1000N1000$mu.medium,
                                   rbind(sd(Scenario7Sim1000N1000$theta[,1]),sd(Scenario7Sim1000N1000$theta[,2])),
                                   Scenario7Sim1000N1000$mu.stdev,AlocScenario7Sim1000N1000$PropAllocation);TableScenario7Sim1000N1000

####=======================================================
#### Scenario 8 - Number of vertices = 2000 - set.seed(13)
####=======================================================
set.seed(13)
Scenario8Sim1000N2000 = EMsimulations(n = 2000, C = 2, nsim = 1000,
                                      theta.real = c(0.3,0.7),
                                      mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                      mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                      mu1.initial = c(0.2,0.2,0.7,0.3,0.1),
                                      mu2.initial = c(0.1,0.5,0.2,0.7,0.7));
AlocScenario8Sim1000N2000 = AlocLikelihoodIndC2(n = 1000, theta.real = c(0.30,0.70),
                                                mu1.real = c(0.00,0.10,0.40,0.10,0.00),
                                                mu2.real = c(0.30,0.40,0.10,0.90,0.80),
                                                mu1.estimate = c(Scenario8Sim1000N2000$mu.medium[1,]),
                                                mu2.estimate = c(Scenario8Sim1000N2000$mu.medium[2,]));
TableScenario8Sim1000N2000 = cbind(rbind(mean(Scenario8Sim1000N2000$theta[,1]),
                                         mean(Scenario8Sim1000N2000$theta[,2])),
                                   Scenario8Sim1000N2000$mu.medium,
                                   rbind(sd(Scenario8Sim1000N2000$theta[,1]),sd(Scenario8Sim1000N2000$theta[,2])),
                                   Scenario8Sim1000N2000$mu.stdev,AlocScenario8Sim1000N2000$PropAllocation);TableScenario8Sim1000N2000

####============================
#### N - Construction of Tables
####============================
Table8_appendix = rbind(TableScenario1Sim1000N20,TableScenario2Sim1000N60,TableScenario3Sim1000N100,
                        TableScenario4Sim1000N150,TableScenario5Sim1000N200,TableScenario6Sim1000N500,
                        TableScenario7Sim1000N1000,TableScenario8Sim1000N2000)

####====================================
#### C - set.seed(13)
#### Changing the number of Communities
####====================================
####=====================================================================
#### Scenario 1 - Number of simulations = 1000 - 3 Communities - N = 100
####=====================================================================
set.seed(13);Scenario1Theta3 = EMsimulations3Communities(n = 100, C = 3, nsim = 1000,
                                                         theta.real = c(0.40,0.20,0.40),
                                                         mu1.real = c(0.00,0.80,0.40), 
                                                         mu2.real = c(0.70,0.30,0.10),
                                                         mu3.real = c(0.20,0.10,0.60),
                                                         mu1.initial = c(0.10,0.60,0.30), 
                                                         mu2.initial = c(0.80,0.15,0.20),
                                                         mu3.initial = c(0.10,0.20,0.50));
AlocScenario1Theta3 = AlocLikelihoodIndC3(n = 100, theta.real = c(0.40,0.20,0.40),
                                          mu1.real = c(0.00,0.80,0.40), 
                                          mu2.real = c(0.70,0.30,0.10),
                                          mu3.real = c(0.20,0.10,0.60),
                                          mu1.estimate = c(Scenario1Theta3$mu.medium[1,]),
                                          mu2.estimate = c(Scenario1Theta3$mu.medium[2,]),
                                          mu3.estimate = c(Scenario1Theta3$mu.medium[3,]));
TableScenario1Theta3 = cbind(rbind(mean(Scenario1Theta3$theta[,1]),
                                   mean(Scenario1Theta3$theta[,2]),
                                   mean(Scenario1Theta3$theta[,3])),
                             Scenario1Theta3$mu.medium,
                             rbind(sd(Scenario1Theta3$theta[,1]),
                                   sd(Scenario1Theta3$theta[,2]),
                                   sd(Scenario1Theta3$theta[,3])),
                             Scenario1Theta3$mu.stdev,AlocScenario1Theta3$PropAllocation);TableScenario1Theta3

####=====================================================================
#### Scenario 2 - Number of simulations = 1000 - 4 Communities - N = 100
####=====================================================================
set.seed(13);
Scenario2Theta4 = EMsimulations4Communities(n = 100, C = 4, nsim = 1000,
                                            theta.real = c(0.25,0.25,0.25,0.25),
                                            mu1.real = c(0.10,0.70,0.10), 
                                            mu2.real = c(0.70,0.20,0.50),
                                            mu3.real = c(0.00,0.20,0.90),
                                            mu4.real = c(0.33,0.33,0.33),
                                            mu1.initial = c(0.20,0.60,0.20), 
                                            mu2.initial = c(0.80,0.30,0.30),
                                            mu3.initial = c(0.10,0.40,0.70),
                                            mu4.initial = c(0.50,0.40,0.30));
AlocScenario2Theta4 = AlocLikelihoodIndC4(n = 100, theta.real = c(0.25,0.25,0.25,0.25),
                                          mu1.real = c(0.10,0.70,0.10), 
                                          mu2.real = c(0.70,0.20,0.50),
                                          mu3.real = c(0.00,0.20,0.90),
                                          mu4.real = c(0.33,0.33,0.33),
                                          mu1.estimate = c(Scenario2Theta4$mu.medium[1,]),
                                          mu2.estimate = c(Scenario2Theta4$mu.medium[2,]),
                                          mu3.estimate = c(Scenario2Theta4$mu.medium[3,]),
                                          mu4.estimate = c(Scenario2Theta4$mu.medium[4,]));
TableScenario2Theta4 = cbind(rbind(mean(Scenario2Theta4$theta[,1]),
                                   mean(Scenario2Theta4$theta[,2]),
                                   mean(Scenario2Theta4$theta[,3]),
                                   mean(Scenario2Theta4$theta[,4])),
                             Scenario2Theta4$mu.medium,
                             rbind(sd(Scenario2Theta4$theta[,1]),
                                   sd(Scenario2Theta4$theta[,2]),
                                   sd(Scenario2Theta4$theta[,3]),
                                   sd(Scenario2Theta4$theta[,4])),
                             Scenario2Theta4$mu.stdev,AlocScenario2Theta4$PropAllocation);TableScenario2Theta4

####=====================================================================
#### Scenario 3 - Number of simulations = 1000 - 5 Communities - N = 100
####=====================================================================
set.seed(13);
Scenario3Theta5 = EMsimulations5Communities(n = 100, C = 5, nsim = 1000,
                                            theta.real = c(0.10,0.20,0.10,0.30,0.30),
                                            mu1.real = c(0.00,0.90,0.50), 
                                            mu2.real = c(0.30,0.00,0.20),
                                            mu3.real = c(0.10,0.05,0.60),
                                            mu4.real = c(0.80,0.20,0.00),
                                            mu5.real = c(0.20,0.00,0.90),
                                            mu1.initial = c(0.10,0.80,0.40), 
                                            mu2.initial = c(0.10,0.10,0.10),
                                            mu3.initial = c(0.20,0.50,0.50),
                                            mu4.initial = c(0.60,0.20,0.20),
                                            mu5.initial = c(0.10,0.20,0.70));
AlocScenario3Theta5 = AlocLikelihoodIndC5(n = 100, theta.real = c(0.10,0.20,0.10,0.30,0.30),
                                          mu1.real = c(0.00,0.90,0.50), 
                                          mu2.real = c(0.30,0.50,0.20),
                                          mu3.real = c(0.10,0.05,0.60),
                                          mu4.real = c(0.80,0.20,0.00),
                                          mu5.real = c(0.20,0.00,0.80),
                                          mu1.estimate = c(Scenario3Theta5$mu.medium[1,]),
                                          mu2.estimate = c(Scenario3Theta5$mu.medium[2,]),
                                          mu3.estimate = c(Scenario3Theta5$mu.medium[3,]),
                                          mu4.estimate = c(Scenario3Theta5$mu.medium[4,]),
                                          mu5.estimate = c(Scenario3Theta5$mu.medium[5,]));
TableScenario3Theta5 = cbind(rbind(mean(Scenario3Theta5$theta[,1]),
                                   mean(Scenario3Theta5$theta[,2]),
                                   mean(Scenario3Theta5$theta[,3]),
                                   mean(Scenario3Theta5$theta[,4]),
                                   mean(Scenario3Theta5$theta[,5])),
                             Scenario3Theta5$mu.medium,
                             rbind(sd(Scenario3Theta5$theta[,1]),
                                   sd(Scenario3Theta5$theta[,2]),
                                   sd(Scenario3Theta5$theta[,3]),
                                   sd(Scenario3Theta5$theta[,4]),
                                   sd(Scenario3Theta5$theta[,5])),
                             Scenario3Theta5$mu.stdev,AlocScenario3Theta5$PropAllocation);TableScenario3Theta5

####============================
#### C - Construction of Tables
####============================
Table9_appendix = rbind(TableScenario1Theta3,TableScenario2Theta4,TableScenario3Theta5)

####==============================
#### m - set.seed(13)
#### Changing the number of words
####==============================
####=================================================
#### Scenario 1 - Number of words = 2 - set.seed(13)
####=================================================
set.seed(13);Scenario1M2 = EMsimulations(n = 100, C = 2, nsim = 1000,
                                         theta.real = c(0.4,0.6),
                                         mu1.real = c(0.70,0.40), 
                                         mu2.real = c(0.00,0.90),
                                         mu1.initial = c(0.80,0.20), 
                                         mu2.initial = c(0.50,0.10));
AlocTableScenario1M2 = AlocLikelihoodIndC2(n = 100, theta.real = c(0.4,0.6),
                                           mu1.real = c(0.70,0.40), 
                                           mu2.real = c(0.00,0.90),
                                           mu1.estimate = c(Scenario1M2$mu.medium[1,]),
                                           mu2.estimate = c(Scenario1M2$mu.medium[2,]));
TableScenario1M2 = cbind(rbind(mean(Scenario1M2$theta[,1]),
                               mean(Scenario1M2$theta[,2])),
                         Scenario1M2$mu.medium,
                         rbind(sd(Scenario1M2$theta[,1]),sd(Scenario1M2$theta[,2])),
                         Scenario1M2$mu.stdev,AlocTableScenario1M2$PropAllocation);TableScenario1M2

####=================================================
#### Scenario 2 - Number of words = 3 - set.seed(13)
####=================================================
set.seed(13);
Scenario2M3 = EMsimulations(n = 100, C = 2, nsim = 1000,
                            theta.real = c(0.8,0.2),
                            mu1.real = c(0.00,0.10,0.50), 
                            mu2.real = c(0.40,0.70,0.05),
                            mu1.initial = c(0.20,0.20,0.20), 
                            mu2.initial = c(0.50,0.90,0.10));
AlocTableScenario2M3 = AlocLikelihoodIndC2(n = 100, theta.real = c(0.8,0.2),
                                           mu1.real = c(0.00,0.10,0.50), 
                                           mu2.real = c(0.40,0.70,0.30),
                                           mu1.estimate = c(Scenario2M3$mu.medium[1,]),
                                           mu2.estimate = c(Scenario2M3$mu.medium[2,]));
TableScenario2M3 = cbind(rbind(mean(Scenario2M3$theta[,1]),
                               mean(Scenario2M3$theta[,2])),
                         Scenario2M3$mu.medium,
                         rbind(sd(Scenario2M3$theta[,1]),sd(Scenario2M3$theta[,2])),
                         Scenario2M3$mu.stdev,AlocTableScenario2M3$PropAllocation);TableScenario2M3

####=================================================
#### Scenario 3 - Number of words = 4 - set.seed(13)
####=================================================
set.seed(13);
Scenario3M4 = EMsimulations(n = 100, C = 2, nsim = 1000,
                            theta.real = c(0.30,0.70),
                            mu1.real = c(0.30,0.60,0.40,0.80), 
                            mu2.real = c(0.90,0.10,0.10,0.00),
                            mu1.initial = c(0.50,0.50,0.50,0.50), 
                            mu2.initial = c(0.30,0.30,0.30,0.30));
AlocTableScenario3M4 = AlocLikelihoodIndC2(n = 100, theta.real = c(0.30,0.70),
                                           mu1.real = c(0.30,0.60,0.40,0.80), 
                                           mu2.real = c(0.90,0.10,0.10,0.00),
                                           mu1.estimate = c(Scenario3M4$mu.medium[1,]),
                                           mu2.estimate = c(Scenario3M4$mu.medium[2,]));
TableScenario3M4 = cbind(rbind(mean(Scenario3M4$theta[,1]),
                               mean(Scenario3M4$theta[,2])),
                         Scenario3M4$mu.medium,
                         rbind(sd(Scenario3M4$theta[,1]),sd(Scenario3M4$theta[,2])),
                         Scenario3M4$mu.stdev,AlocTableScenario3M4$PropAllocation);TableScenario3M4

####=================================================
#### Scenario 4 - Number of words = 6 - set.seed(13)
####=================================================
set.seed(13);
Scenario4M6 = EMsimulations(n = 100, C = 2, nsim = 1000,
                            theta.real = c(0.90,0.10),
                            mu1.real = c(0.10,0.20,0.30,0.40,0.50,0.60), 
                            mu2.real = c(0.80,0.70,0.60,0.00,0.00,0.10),
                            mu1.initial = c(0.20,0.30,0.50,0.30,0.70,0.90), 
                            mu2.initial = c(0.50,0.50,0.50,0.10,0.10,0.10));
AlocTableScenario4M6 = AlocLikelihoodIndC2(n = 100, theta.real = c(0.90,0.10),
                                           mu1.real = c(0.10,0.20,0.30,0.40,0.50,0.60), 
                                           mu2.real = c(0.80,0.70,0.60,0.00,0.00,0.10),
                                           mu1.estimate = c(Scenario4M6$mu.medium[1,]),
                                           mu2.estimate = c(Scenario4M6$mu.medium[2,]));
TableScenario4M6 = cbind(rbind(mean(Scenario4M6$theta[,1]),
                               mean(Scenario4M6$theta[,2])),
                         Scenario4M6$mu.medium,
                         rbind(sd(Scenario4M6$theta[,1]),sd(Scenario4M6$theta[,2])),
                         Scenario4M6$mu.stdev,AlocTableScenario4M6$PropAllocation);TableScenario4M6


####=================================================
#### Scenario 5 - Number of words = 7 - set.seed(13)
####=================================================
set.seed(13);
Scenario5M7 = EMsimulations(n = 100, C = 2, nsim = 1000,
                            theta.real = c(0.10,0.90),
                            mu1.real = c(0.10,0.10,0.10,0.10,0.10,0.10,0.10), 
                            mu2.real = c(0.50,0.20,0.10,0.70,1.00,0.40,0.60),
                            mu1.initial = c(0.30,0.40,0.20,0.10,0.30,0.10,0.50), 
                            mu2.initial = c(0.20,0.05,0.40,0.50,0.10,0.50,0.80));
AlocTableScenario5M7 = AlocLikelihoodIndC2(n = 100, theta.real = c(0.10,0.90),
                                           mu1.real = c(0.10,0.10,0.10,0.10,0.10,0.10,0.10), 
                                           mu2.real = c(0.50,0.20,0.10,0.70,1.00,0.40,0.60),
                                           mu1.estimate = c(Scenario5M7$mu.medium[1,]),
                                           mu2.estimate = c(Scenario5M7$mu.medium[2,]));
TableScenario5M7 = cbind(rbind(mean(Scenario5M7$theta[,1]),
                               mean(Scenario5M7$theta[,2])),
                         Scenario5M7$mu.medium,
                         rbind(sd(Scenario5M7$theta[,1]),sd(Scenario5M7$theta[,2])),
                         Scenario5M7$mu.stdev,AlocTableScenario5M7$PropAllocation);TableScenario5M7

####=================================================
#### Scenario 6 - Number of words = 8 - set.seed(13)
####=================================================
set.seed(13);
Scenario6M8 = EMsimulations(n = 100, C = 2, nsim = 1000,
                            theta.real = c(0.50,0.50),
                            mu1.real = c(0.60,0.00,0.20,0.50,0.30,0.00,0.70,0.10), 
                            mu2.real = c(0.10,0.90,0.50,0.10,0.90,0.40,0.00,0.60),
                            mu1.initial = c(0.20,0.10,0.30,0.60,0.10,0.20,0.40,0.20), 
                            mu2.initial = c(0.30,0.80,0.30,0.20,0.40,0.10,0.05,0.30));
AlocTableScenario6M8 = AlocLikelihoodIndC2(n = 100, theta.real = c(0.50,0.50),
                                           mu1.real = c(0.60,0.00,0.20,0.50,0.30,0.00,0.70,0.10), 
                                           mu2.real = c(0.10,0.90,0.50,0.10,0.90,0.40,0.00,0.60),
                                           mu1.estimate = c(Scenario6M8$mu.medium[1,]),
                                           mu2.estimate = c(Scenario6M8$mu.medium[2,]));
TableScenario6M8 = cbind(rbind(mean(Scenario6M8$theta[,1]),
                               mean(Scenario6M8$theta[,2])),
                         Scenario6M8$mu.medium,
                         rbind(sd(Scenario6M8$theta[,1]),sd(Scenario6M8$theta[,2])),
                         Scenario6M8$mu.stdev,AlocTableScenario6M8$PropAllocation);TableScenario6M8

####=================================================
#### Scenario 7 - Number of words = 9 - set.seed(13)
####=================================================
set.seed(13);
Scenario7M9 = EMsimulations(n = 100, C = 2, nsim = 1000,
                            theta.real = c(0.30,0.70),
                            mu1.real = c(0.20,0.00,0.70,0.50,0.30,0.80,0.05,0.20,0.30), 
                            mu2.real = c(0.60,0.40,0.10,0.00,0.30,0.10,0.90,0.60,0.70),
                            mu1.initial = c(0.10,0.05,0.30,0.40,0.10,0.50,0.10,0.30,0.40),  
                            mu2.initial = c(0.70,0.20,0.20,0.10,0.20,0.20,0.50,0.40,0.50));
AlocTableScenario7M9 = AlocLikelihoodIndC2(n = 100, theta.real = c(0.30,0.70),
                                           mu1.real = c(0.20,0.00,0.70,0.50,0.30,0.80,0.05,0.20,0.30), 
                                           mu2.real = c(0.60,0.40,0.10,0.00,0.30,0.10,0.90,0.60,0.70),
                                           mu1.estimate = c(Scenario7M9$mu.medium[1,]),
                                           mu2.estimate = c(Scenario7M9$mu.medium[2,]));
TableScenario7M9 = cbind(rbind(mean(Scenario7M9$theta[,1]),
                               mean(Scenario7M9$theta[,2])),
                         Scenario7M9$mu.medium,
                         rbind(sd(Scenario7M9$theta[,1]),sd(Scenario7M9$theta[,2])),
                         Scenario7M9$mu.stdev,AlocTableScenario7M9$PropAllocation);TableScenario7M9

####============================
#### m - Construction of Tables
####============================
Table10.1_appendix = TableScenario1M2
Table10.2_appendix = TableScenario2M3
Table10.3_appendix = TableScenario3M4
Table10.4_appendix = TableScenario4M6
Table10.5_appendix = TableScenario5M7
Table10.6_appendix = TableScenario6M8
Table10.7_appendix = TableScenario7M9
Table10.8_appendix = TableScenario8M10
