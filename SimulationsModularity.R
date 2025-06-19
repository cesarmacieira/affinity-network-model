####========================
#### Affinity Network Model
#### César Macieira
####========================
rm(list=ls(all=T))
set.seed(13)

####===========
#### Libraries
####===========
if(!require(openxlsx)){ install.packages("openxlsx"); require(openxlsx)}
if(!require(igraph)){ install.packages("igraph"); require(igraph)}
if(!require(tidyverse)){ install.packages("tidyverse"); require(tidyverse)}
if(!require(tm)){ install.packages("tm", dependencies = T); require(tm) }
if(!require(SnowballC)){ install.packages("SnowballC"); require(SnowballC) }
if(!require(gtools)){ install.packages("gtools"); require(gtools) }

####===========
#### Functions
####===========
AlgorithmEM = function(U_Matrix,C,tol = 1e-8,mu.initial_guess){
  m = ncol(U_Matrix)
  n = nrow(U_Matrix)
  mu.estimate = mu.initial_guess
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
      p1[[k]] = (mumtx[[k]]^U_Matrix) * ((1-mumtx[[k]])^(1-U_Matrix))
      for(i in 1:n){
        p2[i,k] = prod(p1[[k]][i,])
      }
    }
    ptheta = matrix(rep(theta.estimate,n),ncol = C,byrow = T)
    p3 = p2*ptheta # Calculate the numerator of Tik
    # Calculate Tik
    Tik = p3/matrix(rep(rowSums(p3),C),ncol = C)
    # Calculate Q(µ|µ(t))
    parte1 = sum(Tik*log(ptheta))
    parte2 = parte3 = 0 
    for(k in 1:C){
      parte2 = parte2 + sum((Tik*rowSums(log(mumtx[[k]]+1e-300)*U_Matrix))[,k], na.rm = T)
      parte3 = parte3 + sum((Tik*rowSums(log(1-mumtx[[k]]+1e-300)*(1-U_Matrix)))[,k], na.rm = T)
    }
    E = parte1 + parte2 + parte3
    # Estimate θ
    theta.estimate = colMeans(Tik)
    # Estimate µ
    Tikmtx = list()
    for(k in 1:C){
      Tikmtx[[k]] = replicate(m,Tik[,k])
    }
    mu.estimate = matrix(nrow = C, ncol = m)
    for(k in 1:C){
      mu.estimate[k,] = t(colSums(Tikmtx[[k]]*U_Matrix))/colSums(Tikmtx[[k]])
    }
    iterations = iterations + 1
  }
  output = list()
  output$iterations = iterations
  output$loglikelihood = E
  output$theta = thetae
  output$mu = round(mue,3)
  #output$fuzzy = Tik
  return(output)
}
AloclikelihoodIndC2 = function(U_Matrix,Thetas,mu1.estimate,mu2.estimate){
  mu.estimate = rbind(mu1.estimate,mu2.estimate)
  MatrixLikelihoods = matrix(NA, ncol = dim(U_Matrix)[1], 
                             nrow = length(Thetas))
  for(k in 1:length(Thetas)){
    for(i in 1:dim(U_Matrix)[1]){
      MatrixLikelihoods[k,i] = prod(mu.estimate[k,]^U_Matrix[i,])
    }
  }
  Community1 = Community2 = 0
  Allocation = matrix(nrow = 1, ncol = dim(U_Matrix)[1]) 
  for(i in 1:dim(U_Matrix)[1]){
    if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i]){
      sorteioempate = sample(x = c(1,2), size = 1, prob = c(Thetas[1],Thetas[2]))
      if(sorteioempate == 1){
        Community1 = Community1 + 1
        Allocation[1,i] = 1
      }else if(sorteioempate == 2){
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
  Results = list()
  Results$MatrixLikelihoods = MatrixLikelihoods
  Results$Allocation = Allocation
  Results$Theta.estimates = Thetas
  Results$PropAllocation = rbind(Community1/dim(U_Matrix)[1],Community2/dim(U_Matrix)[1])
  return(Results)
}
AlocLikelihoodIndC3 = function(U_Matrix,Thetas,mu1.estimate,mu2.estimate,mu3.estimate){
  mu.estimate = rbind(mu1.estimate,mu2.estimate,mu3.estimate)
  MatrixLikelihoods = matrix(NA, ncol = dim(U_Matrix)[1], nrow = length(Thetas))
  for(k in 1:length(Thetas)){
    for(i in 1:dim(U_Matrix)[1]){
      MatrixLikelihoods[k,i] = prod(mu.estimate[k,]^U_Matrix[i,])
    }
  }
  Community1 = Community2 = Community3 = 0
  Allocation = matrix(nrow = 1, ncol = dim(U_Matrix)[1]) 
  for(i in 1:dim(U_Matrix)[1]){
    if(MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i] & MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i]){
      choicedraw = sample(x = c(1,2,3), size = 1, prob = c(Thetas[1],Thetas[2],Thetas[3]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[2,i]){
      choicedraw = sample(x = c(1,2), size = 1, prob = c(Thetas[1],Thetas[2]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1  
      }
      else if(choicedraw == 2){ 
        Community2 = Community2 + 1; Allocation[1,i] = 2 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 1 & MatrixLikelihoods[1,i] == MatrixLikelihoods[3,i]){
      choicedraw = sample(x = c(1,3), size = 1, prob = c(Thetas[1],Thetas[3]))
      if(choicedraw == 1){ 
        Community1 = Community1 + 1; Allocation[1,i] = 1 
      }
      else if(choicedraw == 3){ 
        Community3 = Community3 + 1; Allocation[1,i] = 3 
      }
    }else if(which.max(MatrixLikelihoods[,i]) == 2 & MatrixLikelihoods[2,i] == MatrixLikelihoods[3,i]){
      choicedraw = sample(x = c(2,3), size = 1, prob = c(Thetas[2],Thetas[3]))
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
  results$PropAllocation = rbind(Community1/dim(U_Matrix)[1],Community2/dim(U_Matrix)[1],Community3/dim(U_Matrix)[1])
  return(results)
  
}
SimulationsC2EM = function(C,nsim,Matriz.real,mu1.initial,mu2.initial){
  U_Matrix = Matriz.real
  theta.sim = matrix(ncol = C, nrow = nsim)
  mu.sim = list()
  likelihood = matrix(ncol = 1, nrow = nsim)
  for(sim in 1:nsim){ # Execute simulations
    mu.initial_guess = rbind(mu1.initial,mu2.initial)
    EstimacaoEM = AlgorithmEM(U_Matrix,C,mu.initial_guess=rbind(mu1.initial,mu2.initial))
    for(k in 1:C){
      theta.sim[sim,k] = EstimacaoEM$theta[k]
      mu.sim[[sim]] = EstimacaoEM$mu
    }
  }
  mu.mean = matrix(nrow = C, ncol = ncol(U_Matrix))
  for(i in 1:C){
    for(j in 1:ncol(U_Matrix)){
      mu.sum = 0
      for(sim in 1:nsim){
        mu.sum = mu.sum + mu.sim[[sim]][i,j]
      }
      mu.mean[i,j] = mu.sum/nsim
    }
  }
  parameters.est = list()
  parameters.est$U_Matrix = U_Matrix
  parameters.est$likelihood = likelihood
  parameters.est$theta = theta.sim
  parameters.est$mu = mu.sim
  parameters.est$mu.mean = mu.mean
  return(parameters.est)
}
SimulationsC3EM = function(C,nsim,Matriz.real,mu1.initial,mu2.initial,mu3.initial){
  U_Matrix = Matriz.real
  theta.sim = matrix(ncol = C, nrow = nsim)
  mu.sim = list()
  likelihood = matrix(ncol = 1, nrow = nsim)
  for(sim in 1:nsim){ # Execute simulations
    mu.initial_guess = rbind(mu1.initial,mu2.initial,mu3.initial)
    EstimacaoEM = AlgorithmEM(U_Matrix,C,mu.initial_guess=rbind(mu1.initial,mu2.initial,mu3.initial))
    for(k in 1:C){
      theta.sim[sim,k] = EstimacaoEM$theta[k]
      mu.sim[[sim]] = EstimacaoEM$mu
    }
  }
  mu.mean = matrix(nrow = C, ncol = ncol(U_Matrix))
  for(i in 1:C){
    for(j in 1:ncol(U_Matrix)){
      mu.sum = 0
      for(sim in 1:nsim){
        mu.sum = mu.sum + mu.sim[[sim]][i,j]
      }
      mu.mean[i,j] = mu.sum/nsim
    }
  }
  parameters.est = list()
  parameters.est$U_Matrix = U_Matrix
  parameters.est$theta = theta.sim
  parameters.est$mu = mu.sim
  parameters.est$mu.mean = mu.mean
  return(parameters.est)
}

####============
#### Modularity
####============
# 1. Updated parameters
n_individuals = 100
n_words = 20
p_ij = matrix(0, nrow = n_individuals, ncol = n_words)
cluster_real = rep(NA, n_individuals)
number_of_people_group_1 = 1:33
number_of_people_group_2 = 34:66
number_of_people_group_3 = 67:100

# 2. Fixed measures by Group
mu1 = c(runif(8, 0.7, 1), runif(12, 0.0, 0.1))
mu2 = c(runif(8, 0.1, 0.2), runif(7, 0.8, 1), runif(5, 0.1, 0.2))
mu3 = c(runif(15, 0.1, 0.3), runif(5, 0.9, 1))

# 3. Group assignment (same number of individuals per Group, adjusted)
for (i in 1:n_individuals) {
  if (i <= number_of_people_group_1[length(number_of_people_group_1)]) {
    p_ij[i, ] = mu1
    cluster_real[i] = 1
  } else if (i <= number_of_people_group_2[length(number_of_people_group_2)]) {
    p_ij[i, ] = mu2
    cluster_real[i] = 2
  } else {
    p_ij[i, ] = mu3
    cluster_real[i] = 3
  }
}

# 4. Generate U_ij
U_ij = matrix(runif(n_individuals * n_words) < p_ij, nrow = n_individuals, ncol = n_words) * 1

# 5. Create bipartite network and design
g_bip = graph.incidence(U_ij)
proj = bipartite_projection(g_bip)$proj1
V(proj)$Group = cluster_real

# 6. Modularity
mod = modularity(proj, membership = cluster_real)
print(mod)

# 7. Detect Communities with Louvain method (can be changed to another)
Communities_mod = cluster_louvain(proj)

# 8. Get Allocation vector (Community detected for each person)
aloc_modularity = Communities_mod$membership

####=====
#### BIC
####=====
# Compute relative frequencies of words for each group
freq_relativa_Group1 = colSums(U_ij[number_of_people_group_1, ]) / sum(U_ij[number_of_people_group_1, ])
freq_relativa_Group2 = colSums(U_ij[number_of_people_group_2, ]) / sum(U_ij[number_of_people_group_2, ])
freq_relativa_Group3 = colSums(U_ij[number_of_people_group_3, ]) / sum(U_ij[number_of_people_group_3, ])

# Number of words (dimensions) per group
kGroup1 = length(freq_relativa_Group1)
kGroup2 = length(freq_relativa_Group2)
kGroup3 = length(freq_relativa_Group3)

# Number of individuals per group
nGroup1 = nrow(U_ij[number_of_people_group_1, ])
nGroup2 = nrow(U_ij[number_of_people_group_2, ])
nGroup3 = nrow(U_ij[number_of_people_group_3, ])

####====================
#### Nº Communities = 1
####====================
# Group 1
theta1.Group1.C1 <- 1  # Initial mixing proportion for 1-component model (fixed at 1)
sumGroup1C1BIC <- matrix(ncol = 1, nrow = nrow(U_ij[number_of_people_group_1, ]))

# Log-likelihood computation for each individual in Group 1
for (i in 1:nrow(U_ij[number_of_people_group_1, ])) {
  sumGroup1C1BIC[i, 1] <- log(sum(theta1.Group1.C1 * (freq_relativa_Group1 ^ U_ij[i, ])))
}

# Number of free parameters = 2 * number of words - 2
deltaC1Group1 <- (2 * kGroup1) - 2

# BIC for Group 1 with 1 component
BICGroup1C1 <- (2 * sum(sumGroup1C1BIC[, 1])) - (deltaC1Group1 * log(nGroup1))
BICGroup1C1

# Group 2
theta1.Group2.C1 <- 1
sumGroup2C1BIC <- matrix(ncol = 1, nrow = nrow(U_ij[number_of_people_group_2, ]))

for (i in 1:nrow(U_ij[number_of_people_group_2, ])) {
  sumGroup2C1BIC[i, 1] <- log(sum(theta1.Group2.C1 * (freq_relativa_Group2 ^ U_ij[i, ])))
}

deltaC1Group2 <- (2 * kGroup2) - 2
BICGroup2C1 <- (2 * sum(sumGroup2C1BIC[, 1])) - (deltaC1Group2 * log(nGroup2))
BICGroup2C1

# Group 3
theta1.Group3.C1 <- 1
sumGroup3C1BIC <- matrix(ncol = 1, nrow = nrow(U_ij[number_of_people_group_3, ]))

for (i in 1:nrow(U_ij[number_of_people_group_3, ])) {
  sumGroup3C1BIC[i, 1] <- log(sum(theta1.Group3.C1 * (freq_relativa_Group3 ^ U_ij[i, ])))
}

deltaC1Group3 <- (2 * kGroup3) - 2
BICGroup3C1 <- (2 * sum(sumGroup3C1BIC[, 1])) - (deltaC1Group3 * log(nGroup3))
BICGroup3C1

####====================
#### Nº Communities = 2
####====================
# Group 1
# Set seed for reproducibility of initial mu1 for C=2 model
set.seed(13)
Group1C2initialC1 = rdirichlet(n = 1, alpha = rep(1, dim(U_ij[number_of_people_group_1, ])[2]))

# Set another seed for reproducibility of initial mu2 for C=2 model
set.seed(1968)
Group1C2initialC2 = rdirichlet(n = 1, alpha = rep(1, dim(U_ij[number_of_people_group_1, ])[2]))

# Dictionary to name word columns
dicionarioGroup1 = c(seq(1:20))

# Run EM algorithm for C = 2 clusters on Group 1 data
Group1C2BIC = SimulationsC2EM(
  C = 2,
  nsim = 1,
  Matriz.real = U_ij[number_of_people_group_1, ],
  mu1.initial = Group1C2initialC1,
  mu2.initial = Group1C2initialC2
)

# Combine resulting theta and mu estimates into a single table
TableGroup1C2BIC = cbind(
  rbind(mean(Group1C2BIC$theta[, 1]), mean(Group1C2BIC$theta[, 2])),
  Group1C2BIC$mu.mean
)
colnames(TableGroup1C2BIC) = c("theta", dicionarioGroup1)
TableGroup1C2BIC

# Initialize matrix to store log-likelihood per individual
sumGroup1C2BIC = matrix(ncol = 1, nrow = dim(U_ij[number_of_people_group_1, ])[1])

# Extract theta and mu estimates from the EM output
theta1.Group1.C2 = TableGroup1C2BIC[1, 1]
theta2.Group1.C2 = TableGroup1C2BIC[2, 1]

mu1.estimated.Group1.C2 = TableGroup1C2BIC[1, 2:ncol(TableGroup1C2BIC)]
mu1.standardized.Group1.C2 = mu1.estimated.Group1.C2 / sum(mu1.estimated.Group1.C2)

mu2.estimated.Group1.C2 = TableGroup1C2BIC[2, 2:ncol(TableGroup1C2BIC)]
mu2.standardized.Group1.C2 = mu2.estimated.Group1.C2 / sum(mu2.estimated.Group1.C2)

# Compute the log-likelihood for each individual using the mixture model
for (i in 1:nrow(U_ij[number_of_people_group_1, ])) {
  sumGroup1C2BIC[i, 1] = log(
    (theta1.Group1.C2 * sum(mu1.standardized.Group1.C2^U_ij[i, ])) +
      (theta2.Group1.C2 * sum(mu2.standardized.Group1.C2^U_ij[i, ]))
  )
}
sumGroup1C2BIC

# Degrees of freedom: 2*k - 1 (theta adds one parameter vs C=1)
deltaC2Group1 = (2 * kGroup1) - 1

# Compute BIC for Group 1 with 2 components
BICGroup1C2 = (2 * sum(sumGroup1C2BIC[, 1])) - (deltaC2Group1 * log(nGroup1))
BICGroup1C2

# Group 2
set.seed(13)
Group2C2initialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_2,])[2]))
set.seed(1968)
Group2C2initialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_2,])[2]))

dicionarioGroup2 = c(seq(1:20))

Group2C2BIC = 
  SimulationsC2EM(C = 2, nsim = 1, Matriz.real = U_ij[number_of_people_group_2,],
                  mu1.initial = Group2C2initialC1,
                  mu2.initial = Group2C2initialC2)
TableGroup2C2BIC = cbind(rbind(mean(Group2C2BIC$theta[,1]),mean(Group2C2BIC$theta[,2])),
                         Group2C2BIC$mu.mean)
colnames(TableGroup2C2BIC) = c("theta",dicionarioGroup2);TableGroup2C2BIC

sumGroup2C2BIC = matrix(ncol = 1, nrow = dim(U_ij[number_of_people_group_2,])[1])

theta1.Group2.C2 = TableGroup2C2BIC[1,1]
theta2.Group2.C2 = TableGroup2C2BIC[2,1]

mu1.estimated.Group2.C2 = TableGroup2C2BIC[1,2:dim(TableGroup2C2BIC)[2]]
mu1.standardized.Group2.C2 = mu1.estimated.Group2.C2/sum(TableGroup2C2BIC[1,2:dim(TableGroup2C2BIC)[2]])
mu2.estimated.Group2.C2 = TableGroup2C2BIC[2,2:dim(TableGroup2C2BIC)[2]]
mu2.standardized.Group2.C2 = mu2.estimated.Group2.C2/sum(TableGroup2C2BIC[2,2:dim(TableGroup2C2BIC)[2]])

for(i in 1:dim(U_ij[number_of_people_group_2,])[1]){
  sumGroup2C2BIC[i,1] = 
    log( 
      (theta1.Group2.C2 * sum(mu1.standardized.Group2.C2^U_ij[i,]))+
        (theta2.Group2.C2 * sum(mu2.standardized.Group2.C2^U_ij[i,]))
    )
}
sumGroup2C2BIC
deltaC2Group2 = (2*kGroup2)-1
BICGroup2C2 = (2*sum(sumGroup2C2BIC[,1]))-(deltaC2Group2*log(nGroup2))
BICGroup2C2

# Group 3
set.seed(13)
Group3C2initialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_3,])[2]))
set.seed(1968)
Group3C2initialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_3,])[2]))

dicionarioGroup3 = c(seq(1:20))

Group3C2BIC = 
  SimulationsC2EM(C = 2, nsim = 1, Matriz.real = U_ij[number_of_people_group_3,],
                  mu1.initial = Group3C2initialC1,
                  mu2.initial = Group3C2initialC2)
TableGroup3C2BIC = cbind(rbind(mean(Group3C2BIC$theta[,1]),mean(Group3C2BIC$theta[,2])),
                         Group3C2BIC$mu.mean)
colnames(TableGroup3C2BIC) = c("theta",dicionarioGroup3);TableGroup3C2BIC

sumGroup3C2BIC = matrix(ncol = 1, nrow = dim(U_ij[number_of_people_group_3,])[1])

theta1.Group3.C2 = TableGroup3C2BIC[1,1]
theta2.Group3.C2 = TableGroup3C2BIC[2,1]

mu1.estimated.Group3.C2 = TableGroup3C2BIC[1,2:dim(TableGroup3C2BIC)[2]]
mu1.standardized.Group3.C2 = mu1.estimated.Group3.C2/sum(TableGroup3C2BIC[1,2:dim(TableGroup3C2BIC)[2]])
mu2.estimated.Group3.C2 = TableGroup3C2BIC[2,2:dim(TableGroup3C2BIC)[2]]
mu2.standardized.Group3.C2 = mu2.estimated.Group3.C2/sum(TableGroup3C2BIC[2,2:dim(TableGroup3C2BIC)[2]])

for(i in 1:dim(U_ij[number_of_people_group_3,])[1]){
  sumGroup3C2BIC[i,1] = 
    log( 
      (theta1.Group3.C2 * sum(mu1.standardized.Group3.C2^U_ij[i,]))+
        (theta2.Group3.C2 * sum(mu2.standardized.Group3.C2^U_ij[i,]))
    )
}
sumGroup3C2BIC
deltaC2Group3 = (2*kGroup3)-1
BICGroup3C2 = (2*sum(sumGroup3C2BIC[,1]))-(deltaC2Group3*log(nGroup3))
BICGroup3C2

####====================
#### Nº Communities = 3
####====================
# Group 1
# Set seeds to generate reproducible initial values for the EM algorithm with 3 clusters
set.seed(13)
Group1C3initialC1 = rdirichlet(n = 1, alpha = rep(1, dim(U_ij[number_of_people_group_1, ])[2]))

set.seed(1968)
Group1C3initialC2 = rdirichlet(n = 1, alpha = rep(1, dim(U_ij[number_of_people_group_1, ])[2]))

set.seed(100)
Group1C3initialC3 = rdirichlet(n = 1, alpha = rep(1, dim(U_ij[number_of_people_group_1, ])[2]))

# Run the EM algorithm for C = 3 on Group 1 data
Group1C3BIC = SimulationsC3EM(
  C = 3,
  nsim = 1,
  Matriz.real = U_ij[number_of_people_group_1, ],
  mu1.initial = Group1C3initialC1,
  mu2.initial = Group1C3initialC2,
  mu3.initial = Group1C3initialC3
)

# Combine the estimated theta and mu values into a result table
TableGroup1C3BIC = cbind(
  rbind(
    mean(Group1C3BIC$theta[, 1]),
    mean(Group1C3BIC$theta[, 2]),
    mean(Group1C3BIC$theta[, 3])
  ),
  Group1C3BIC$mu.mean
)
colnames(TableGroup1C3BIC) = c("theta", dicionarioGroup1)
TableGroup1C3BIC

# Initialize log-likelihood matrix for each individual in Group 1
sumGroup1C3BIC = matrix(ncol = 1, nrow = nrow(U_ij[number_of_people_group_1, ]))

# Extract and normalize theta and mu estimates
theta1.Group1.C3 = TableGroup1C3BIC[1, 1]
theta2.Group1.C3 = TableGroup1C3BIC[2, 1]
theta3.Group1.C3 = TableGroup1C3BIC[3, 1]

mu1.estimated.Group1.C3 = TableGroup1C3BIC[1, 2:ncol(TableGroup1C3BIC)]
mu1.standardized.Group1.C3 = mu1.estimated.Group1.C3 / sum(mu1.estimated.Group1.C3)

mu2.estimated.Group1.C3 = TableGroup1C3BIC[2, 2:ncol(TableGroup1C3BIC)]
mu2.standardized.Group1.C3 = mu2.estimated.Group1.C3 / sum(mu2.estimated.Group1.C3)

mu3.estimated.Group1.C3 = TableGroup1C3BIC[3, 2:ncol(TableGroup1C3BIC)]
mu3.standardized.Group1.C3 = mu3.estimated.Group1.C3 / sum(mu3.estimated.Group1.C3)

# Calculate log-likelihood for each individual under the 3-component model
for (i in 1:nrow(U_ij[number_of_people_group_1, ])) {
  sumGroup1C3BIC[i, 1] = log(
    (theta1.Group1.C3 * sum(mu1.standardized.Group1.C3 ^ U_ij[i, ])) +
      (theta2.Group1.C3 * sum(mu2.standardized.Group1.C3 ^ U_ij[i, ])) +
      (theta3.Group1.C3 * sum(mu3.standardized.Group1.C3 ^ U_ij[i, ]))
  )
}
sumGroup1C3BIC

# Degrees of freedom: for C=3, we estimate 3*k parameters for mu and 2 for theta (since they sum to 1), so total = 2*k
deltaC3Group1 = 2 * kGroup1

# Compute BIC for Group 1 under the 3-component model
BICGroup1C3 = (2 * sum(sumGroup1C3BIC[, 1])) - (deltaC3Group1 * log(nGroup1))
BICGroup1C3

# Group 2
set.seed(13)
Group2C3initialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_2,])[2]))
set.seed(1968)
Group2C3initialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_2,])[2]))
set.seed(100)
Group2C3initialC3 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_2,])[2]))

Group2C3BIC = 
  SimulationsC3EM(C = 3, nsim = 1, Matriz.real = U_ij[number_of_people_group_2,],
                  mu1.initial = Group2C3initialC1,
                  mu2.initial = Group2C3initialC2,
                  mu3.initial = Group2C3initialC3)
TableGroup2C3BIC = cbind(rbind(mean(Group2C3BIC$theta[,1]),
                               mean(Group2C3BIC$theta[,2]),
                               mean(Group2C3BIC$theta[,3])),
                         Group2C3BIC$mu.mean)
colnames(TableGroup2C3BIC) = c("theta",dicionarioGroup2);TableGroup2C3BIC

sumGroup2C3BIC = matrix(ncol = 1, nrow = dim(U_ij[number_of_people_group_2,])[1])

theta1.Group2.C3 = TableGroup2C3BIC[1,1]
theta2.Group2.C3 = TableGroup2C3BIC[2,1]
theta3.Group2.C3 = TableGroup2C3BIC[3,1]

mu1.estimated.Group2.C3 = TableGroup2C3BIC[1,2:dim(TableGroup2C3BIC)[2]]
mu1.standardized.Group2.C3 = mu1.estimated.Group2.C3/sum(TableGroup2C3BIC[1,2:dim(TableGroup2C3BIC)[2]])
mu2.estimated.Group2.C3 = TableGroup2C3BIC[2,2:dim(TableGroup2C3BIC)[2]]
mu2.standardized.Group2.C3 = mu2.estimated.Group2.C3/sum(TableGroup2C3BIC[2,2:dim(TableGroup2C3BIC)[2]])
mu3.estimated.Group2.C3 = TableGroup2C3BIC[3,2:dim(TableGroup2C3BIC)[2]]
mu3.standardized.Group2.C3 = mu3.estimated.Group2.C3/sum(TableGroup2C3BIC[3,2:dim(TableGroup2C3BIC)[2]])

for(i in 1:dim(U_ij[number_of_people_group_2,])[1]){
  sumGroup2C3BIC[i,1] = 
    log(
      (theta1.Group2.C3 * sum((mu1.standardized.Group2.C3^U_ij[i,]) ))+
        (theta2.Group2.C3 * sum((mu2.standardized.Group2.C3^U_ij[i,]) ))+
        (theta3.Group2.C3 * sum((mu3.standardized.Group2.C3^U_ij[i,]) ))
    )
}
sumGroup2C3BIC
deltaC3Group2 = (2*kGroup2)
BICGroup2C3 = (2*sum(sumGroup2C3BIC[,1]))-(deltaC3Group2*log(nGroup2))
BICGroup2C3

# Group 3
set.seed(13)
Group3C3initialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_3,])[2]))
set.seed(1968)
Group3C3initialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_3,])[2]))
set.seed(100)
Group3C3initialC3 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[number_of_people_group_3,])[2]))

Group3C3BIC = 
  SimulationsC3EM(C = 3, nsim = 1, Matriz.real = U_ij[number_of_people_group_3,],
                  mu1.initial = Group3C3initialC1,
                  mu2.initial = Group3C3initialC2,
                  mu3.initial = Group3C3initialC3)
TableGroup3C3BIC = cbind(rbind(mean(Group3C3BIC$theta[,1]),
                               mean(Group3C3BIC$theta[,2]),
                               mean(Group3C3BIC$theta[,3])),
                         Group3C3BIC$mu.mean)
colnames(TableGroup3C3BIC) = c("theta",dicionarioGroup3);TableGroup3C3BIC

sumGroup3C3BIC = matrix(ncol = 1, nrow = dim(U_ij[number_of_people_group_3,])[1])

theta1.Group3.C3 = TableGroup3C3BIC[1,1]
theta2.Group3.C3 = TableGroup3C3BIC[2,1]
theta3.Group3.C3 = TableGroup3C3BIC[3,1]

mu1.estimated.Group3.C3 = TableGroup3C3BIC[1,2:dim(TableGroup3C3BIC)[2]]
mu1.standardized.Group3.C3 = mu1.estimated.Group3.C3/sum(TableGroup3C3BIC[1,2:dim(TableGroup3C3BIC)[2]])
mu2.estimated.Group3.C3 = TableGroup3C3BIC[2,2:dim(TableGroup3C3BIC)[2]]
mu2.standardized.Group3.C3 = mu2.estimated.Group3.C3/sum(TableGroup3C3BIC[2,2:dim(TableGroup3C3BIC)[2]])
mu3.estimated.Group3.C3 = TableGroup3C3BIC[3,2:dim(TableGroup3C3BIC)[2]]
mu3.standardized.Group3.C3 = mu3.estimated.Group3.C3/sum(TableGroup3C3BIC[3,2:dim(TableGroup3C3BIC)[2]])

for(i in 1:dim(U_ij[number_of_people_group_3,])[1]){
  sumGroup3C3BIC[i,1] = 
    log(
      (theta1.Group3.C3 * sum((mu1.standardized.Group3.C3^U_ij[i,]) ))+
        (theta2.Group3.C3 * sum((mu2.standardized.Group3.C3^U_ij[i,]) ))+
        (theta3.Group3.C3 * sum((mu3.standardized.Group3.C3^U_ij[i,]) ))
    )
}
sumGroup3C3BIC
deltaC3Group3 = (2*kGroup3)
BICGroup3C3 = (2*sum(sumGroup3C3BIC[,1]))-(deltaC3Group3*log(nGroup3))
BICGroup3C3

####=============================
#### Comparison of results - BIC
####=============================
rbind(BICGroup1C1,BICGroup1C2,BICGroup1C3)
rbind(BICGroup2C1,BICGroup2C2,BICGroup2C3)
rbind(BICGroup3C1,BICGroup3C2,BICGroup3C3)

####=======================
#### Comparison of results
####=======================
# Real group
true_group = c(rep(1, 33), rep(2, 33), rep(3, 34))

# Modularity
aloc_modularity
table(aloc_modularity)

# Likelihood
Results_aloc_EM = AlocLikelihoodIndC3(U_Matrix = U_ij, Thetas = c(0.33,0.33,0.34),
                                      mu1.estimate = mu1, mu2.estimate = mu2, mu3.estimate = mu3)
Results_aloc_EM

df_comparison = data.frame(
  True_Group = true_group,
  Modularity_Cluster = aloc_modularity,
  EM_Cluster = as.vector(Results_aloc_EM$Allocation))

# Direct comparison (assuming label correspondence is correct)
modularity_accuracy <- mean(df_comparison$True_Group == df_comparison$Modularity_Cluster)
em_accuracy <- mean(df_comparison$True_Group == df_comparison$EM_Cluster)

# Print accuracies
cat("Modularity Accuracy:", round(modularity_accuracy * 100, 2), "%\n")
cat("EM Likelihood Accuracy:", round(em_accuracy * 100, 2), "%\n")