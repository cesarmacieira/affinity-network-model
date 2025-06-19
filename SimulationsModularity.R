set.seed(13)
library(igraph)

AlocVerossimilhancaIndC2 <- function(MatrizU,Thetas,
                                     mu1.estimado,mu2.estimado){
  mu.estimado <- rbind(mu1.estimado,mu2.estimado)
  MatrizVerossimilhancas <- matrix(NA, ncol = dim(MatrizU)[1], 
                                   nrow = length(Thetas))
  for(k in 1:length(Thetas)){
    for(i in 1:dim(MatrizU)[1]){
      MatrizVerossimilhancas[k,i] <- prod(mu.estimado[k,]^MatrizU[i,])
    }
  }
  Comunidade1 <- Comunidade2 <- 0
  Alocacao <- matrix(nrow = 1, ncol = dim(MatrizU)[1]) 
  for(i in 1:dim(MatrizU)[1]){
    if(MatrizVerossimilhancas[1,i] == MatrizVerossimilhancas[2,i]){
      sorteioempate <- sample(x = c(1,2), size = 1, prob = c(Thetas[1],Thetas[2]))
      if(sorteioempate == 1){
        Comunidade1 <- Comunidade1 + 1
        Alocacao[1,i] <- 1
      }else if(sorteioempate == 2){
        Comunidade2 <- Comunidade2 + 1
        Alocacao[1,i] <- 2
      }
    }else 
      if(which.max(MatrizVerossimilhancas[,i]) == 1){
        Comunidade1 <- Comunidade1 + 1
        Alocacao[1,i] <- 1
      }else if(which.max(MatrizVerossimilhancas[,i]) == 2){
        Comunidade2 <- Comunidade2 + 1
        Alocacao[1,i] <- 2
      }
  }
  resultados <- list()
  resultados$MatrizVerossimilhancas <- MatrizVerossimilhancas
  resultados$Alocação <- Alocacao
  resultados$Theta.estimados <- Thetas
  resultados$PropAlocação <- rbind(Comunidade1/dim(MatrizU)[1],Comunidade2/dim(MatrizU)[1])
  return(resultados)
}
AlocLikelihoodIndC3 = function(MatrizU,Thetas,
                               mu1.estimate,mu2.estimate,mu3.estimate){
  mu.estimado = rbind(mu1.estimate,mu2.estimate,mu3.estimate)
  MatrixLikelihoods = matrix(NA, ncol = dim(MatrizU)[1], nrow = length(Thetas))
  for(k in 1:length(Thetas)){
    for(i in 1:dim(MatrizU)[1]){
      MatrixLikelihoods[k,i] = prod(mu.estimado[k,]^MatrizU[i,])
    }
  }
  Community1 = Community2 = Community3 = 0
  Allocation = matrix(nrow = 1, ncol = dim(MatrizU)[1]) 
  for(i in 1:dim(MatrizU)[1]){
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
  results$PropAllocation = rbind(Community1/dim(MatrizU)[1],Community2/dim(MatrizU)[1],Community3/dim(MatrizU)[1])
  return(results)
  
}

# 1. Parâmetros atualizados
n_pessoas <- 100
n_palavras <- 20
p_ij <- matrix(0, nrow = n_pessoas, ncol = n_palavras)
cluster_real <- rep(NA, n_pessoas)

# 2. Medidas fixas por grupo
medida1 <- c(runif(8, 0.7, 1), runif(12, 0.0, 0.1))
medida2 <- c(runif(8, 0.1, 0.2), runif(7, 0.8, 1), runif(5, 0.1, 0.2))
medida3 <- c(runif(15, 0.1, 0.3), runif(5, 0.9, 1))

# 3. Atribuição dos grupos (mesmo número de pessoas por grupo, ajustado)
for (i in 1:n_pessoas) {
  if (i <= 33) {
    p_ij[i, ] <- medida1
    cluster_real[i] <- 1
  } else if (i <= 66) {
    p_ij[i, ] <- medida2
    cluster_real[i] <- 2
  } else {
    p_ij[i, ] <- medida3
    cluster_real[i] <- 3
  }
}

# 4. Gerar U_ij
U_ij <- matrix(runif(n_pessoas * n_palavras) < p_ij, nrow = n_pessoas, ncol = n_palavras) * 1

# 5. Criar rede bipartida e projetar
g_bip <- graph.incidence(U_ij)
proj <- bipartite.projection(g_bip)$proj1
V(proj)$grupo <- cluster_real

# 6. Modularidade
mod <- modularity(proj, membership = cluster_real)
print(mod)

# 7. Detectar comunidades com método de Louvain (pode trocar por outro)
comunidades_mod <- cluster_louvain(proj)

# 8. Obter vetor de alocação (comunidade detectada para cada pessoa)
aloc_modularidade <- comunidades_mod$membership

# 9. Ver alocação
print(aloc_modularidade)
table(aloc_modularidade)
table(resultados_aloc_EM$Allocation)

####=====
#### BIC
####=====
freq_relativa_Grupo1 = colSums(U_ij[1:33,]) / sum(U_ij[1:33,])
freq_relativa_Grupo2 = colSums(U_ij[34:66,]) / sum(U_ij[34:66,])
freq_relativa_Grupo3 = colSums(U_ij[67:100,]) / sum(U_ij[67:100,])

kGrupo1 = length(freq_relativa_Grupo1)
nGrupo1 = dim(U_ij[1:33,])[1]

kGrupo2 = length(freq_relativa_Grupo2)
nGrupo2 = dim(U_ij[34:66,])[1]

kGrupo3 = length(freq_relativa_Grupo3)
nGrupo3 = dim(U_ij[67:100,])[1]

####=======================
#### Nº de comunidades = 1
####=======================
# Grupo 1
theta1.Grupo1.C1 = 1
somaGrupo1C1BIC = matrix(ncol = 1, nrow = dim(U_ij[1:33,])[1])
for(i in 1:dim(U_ij[1:33,])[1]){
  somaGrupo1C1BIC[i,1] = log(sum(theta1.Grupo1.C1*(freq_relativa_Grupo1^U_ij[i,])))
}
somaGrupo1C1BIC

deltaC1Grupo1 = (2*kGrupo1)-2
BICGrupo1C1 = (2*sum(somaGrupo1C1BIC[,1]))-(deltaC1Grupo1*log(nGrupo1))
BICGrupo1C1

# Grupo 2
theta1.Grupo2.C1 = 1
somaGrupo2C1BIC = matrix(ncol = 1, nrow = dim(U_ij[34:66,])[1])
for(i in 1:dim(U_ij[34:66,])[1]){
  somaGrupo2C1BIC[i,1] = log(sum(theta1.Grupo2.C1*(freq_relativa_Grupo2^U_ij[i,])))
}
somaGrupo2C1BIC

deltaC1Grupo2 = (2*kGrupo2)-2
BICGrupo2C1 = (2*sum(somaGrupo2C1BIC[,1]))-(deltaC1Grupo2*log(nGrupo2))
BICGrupo2C1

# Grupo 3
theta1.Grupo3.C1 = 1
somaGrupo3C1BIC = matrix(ncol = 1, nrow = dim(U_ij[67:100,])[1])
for(i in 1:dim(U_ij[67:100,])[1]){
  somaGrupo3C1BIC[i,1] = log(sum(theta1.Grupo3.C1*(freq_relativa_Grupo3^U_ij[i,])))
}
somaGrupo3C1BIC

deltaC1Grupo3 = (2*kGrupo3)-2
BICGrupo3C1 = (2*sum(somaGrupo3C1BIC[,1]))-(deltaC1Grupo3*log(nGrupo3))
BICGrupo3C1

####=======================
#### Nº de comunidades = 2
####=======================
# Grupo 1
set.seed(13)
Grupo1C2inicialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[1:33,])[2]))
set.seed(1968)
Grupo1C2inicialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[1:33,])[2]))

dicionarioGrupo1 = c(seq(1:20))

Grupo1C2BIC = 
  SimulacoesC2EM(C = 2, nsim = 1, Matriz.real = U_ij[1:33,],
                 mu1.inicial = Grupo1C2inicialC1,
                 mu2.inicial = Grupo1C2inicialC2)
TabelaGrupo1C2BIC = cbind(rbind(mean(Grupo1C2BIC$theta[,1]),mean(Grupo1C2BIC$theta[,2])),
                           Grupo1C2BIC$mu.media)
colnames(TabelaGrupo1C2BIC) = c("theta",dicionarioGrupo1);TabelaGrupo1C2BIC

somaGrupo1C2BIC = matrix(ncol = 1, nrow = dim(U_ij[1:33,])[1])

theta1.Grupo1.C2 = TabelaGrupo1C2BIC[1,1]
theta2.Grupo1.C2 = TabelaGrupo1C2BIC[2,1]

mu1.estimada.Grupo1.C2 = TabelaGrupo1C2BIC[1,2:dim(TabelaGrupo1C2BIC)[2]]
mu1.padronizada.Grupo1.C2 = mu1.estimada.Grupo1.C2/sum(TabelaGrupo1C2BIC[1,2:dim(TabelaGrupo1C2BIC)[2]])
mu2.estimada.Grupo1.C2 = TabelaGrupo1C2BIC[2,2:dim(TabelaGrupo1C2BIC)[2]]
mu2.padronizada.Grupo1.C2 = mu2.estimada.Grupo1.C2/sum(TabelaGrupo1C2BIC[2,2:dim(TabelaGrupo1C2BIC)[2]])

for(i in 1:dim(U_ij[1:33,])[1]){
  somaGrupo1C2BIC[i,1] = 
    log( 
      (theta1.Grupo1.C2 * sum(mu1.padronizada.Grupo1.C2^U_ij[i,]))+
        (theta2.Grupo1.C2 * sum(mu2.padronizada.Grupo1.C2^U_ij[i,]))
    )
}
somaGrupo1C2BIC
deltaC2Grupo1 = (2*kGrupo1)-1
BICGrupo1C2 = (2*sum(somaGrupo1C2BIC[,1]))-(deltaC2Grupo1*log(nGrupo1))
BICGrupo1C2

# Grupo 2
set.seed(13)
Grupo2C2inicialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[34:66,])[2]))
set.seed(1968)
Grupo2C2inicialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[34:66,])[2]))

dicionarioGrupo2 = c(seq(1:20))

Grupo2C2BIC = 
  SimulacoesC2EM(C = 2, nsim = 1, Matriz.real = U_ij[34:66,],
                 mu1.inicial = Grupo2C2inicialC1,
                 mu2.inicial = Grupo2C2inicialC2)
TabelaGrupo2C2BIC = cbind(rbind(mean(Grupo2C2BIC$theta[,1]),mean(Grupo2C2BIC$theta[,2])),
                           Grupo2C2BIC$mu.media)
colnames(TabelaGrupo2C2BIC) = c("theta",dicionarioGrupo2);TabelaGrupo2C2BIC

somaGrupo2C2BIC = matrix(ncol = 1, nrow = dim(U_ij[34:66,])[1])

theta1.Grupo2.C2 = TabelaGrupo2C2BIC[1,1]
theta2.Grupo2.C2 = TabelaGrupo2C2BIC[2,1]

mu1.estimada.Grupo2.C2 = TabelaGrupo2C2BIC[1,2:dim(TabelaGrupo2C2BIC)[2]]
mu1.padronizada.Grupo2.C2 = mu1.estimada.Grupo2.C2/sum(TabelaGrupo2C2BIC[1,2:dim(TabelaGrupo2C2BIC)[2]])
mu2.estimada.Grupo2.C2 = TabelaGrupo2C2BIC[2,2:dim(TabelaGrupo2C2BIC)[2]]
mu2.padronizada.Grupo2.C2 = mu2.estimada.Grupo2.C2/sum(TabelaGrupo2C2BIC[2,2:dim(TabelaGrupo2C2BIC)[2]])

for(i in 1:dim(U_ij[34:66,])[1]){
  somaGrupo2C2BIC[i,1] = 
    log( 
      (theta1.Grupo2.C2 * sum(mu1.padronizada.Grupo2.C2^U_ij[i,]))+
        (theta2.Grupo2.C2 * sum(mu2.padronizada.Grupo2.C2^U_ij[i,]))
    )
}
somaGrupo2C2BIC
deltaC2Grupo2 = (2*kGrupo2)-1
BICGrupo2C2 = (2*sum(somaGrupo2C2BIC[,1]))-(deltaC2Grupo2*log(nGrupo2))
BICGrupo2C2

# Grupo 3
set.seed(13)
Grupo3C2inicialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[67:100,])[2]))
set.seed(1968)
Grupo3C2inicialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[67:100,])[2]))

dicionarioGrupo3 = c(seq(1:20))

Grupo3C2BIC = 
  SimulacoesC2EM(C = 2, nsim = 1, Matriz.real = U_ij[67:100,],
                 mu1.inicial = Grupo3C2inicialC1,
                 mu2.inicial = Grupo3C2inicialC2)
TabelaGrupo3C2BIC = cbind(rbind(mean(Grupo3C2BIC$theta[,1]),mean(Grupo3C2BIC$theta[,2])),
                           Grupo3C2BIC$mu.media)
colnames(TabelaGrupo3C2BIC) = c("theta",dicionarioGrupo3);TabelaGrupo3C2BIC

somaGrupo3C2BIC = matrix(ncol = 1, nrow = dim(U_ij[67:100,])[1])

theta1.Grupo3.C2 = TabelaGrupo3C2BIC[1,1]
theta2.Grupo3.C2 = TabelaGrupo3C2BIC[2,1]

mu1.estimada.Grupo3.C2 = TabelaGrupo3C2BIC[1,2:dim(TabelaGrupo3C2BIC)[2]]
mu1.padronizada.Grupo3.C2 = mu1.estimada.Grupo3.C2/sum(TabelaGrupo3C2BIC[1,2:dim(TabelaGrupo3C2BIC)[2]])
mu2.estimada.Grupo3.C2 = TabelaGrupo3C2BIC[2,2:dim(TabelaGrupo3C2BIC)[2]]
mu2.padronizada.Grupo3.C2 = mu2.estimada.Grupo3.C2/sum(TabelaGrupo3C2BIC[2,2:dim(TabelaGrupo3C2BIC)[2]])

for(i in 1:dim(U_ij[67:100,])[1]){
  somaGrupo3C2BIC[i,1] = 
    log( 
      (theta1.Grupo3.C2 * sum(mu1.padronizada.Grupo3.C2^U_ij[i,]))+
        (theta2.Grupo3.C2 * sum(mu2.padronizada.Grupo3.C2^U_ij[i,]))
    )
}
somaGrupo3C2BIC
deltaC2Grupo3 = (2*kGrupo3)-1
BICGrupo3C2 = (2*sum(somaGrupo3C2BIC[,1]))-(deltaC2Grupo3*log(nGrupo3))
BICGrupo3C2

####=======================
#### Nº de comunidades = 3
####=======================
# Grupo 1
set.seed(13)
Grupo1C3inicialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[1:33,])[2]))
set.seed(1968)
Grupo1C3inicialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[1:33,])[2]))
set.seed(100)
Grupo1C3inicialC3 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[1:33,])[2]))

#Tempo grande para execução
Grupo1C3BIC = 
  SimulacoesC3EM(C = 3, nsim = 1, Matriz.real = U_ij[1:33,],
                 mu1.inicial = Grupo1C3inicialC1,
                 mu2.inicial = Grupo1C3inicialC2,
                 mu3.inicial = Grupo1C3inicialC3)
TabelaGrupo1C3BIC = cbind(rbind(mean(Grupo1C3BIC$theta[,1]),
                                mean(Grupo1C3BIC$theta[,2]),
                                mean(Grupo1C3BIC$theta[,3])),
                          Grupo1C3BIC$mu.media)
colnames(TabelaGrupo1C3BIC) = c("theta",dicionarioGrupo1);TabelaGrupo1C3BIC

somaGrupo1C3BIC = matrix(ncol = 1, nrow = dim(U_ij[1:33,])[1])

theta1.Grupo1.C3 = TabelaGrupo1C3BIC[1,1]
theta2.Grupo1.C3 = TabelaGrupo1C3BIC[2,1]
theta3.Grupo1.C3 = TabelaGrupo1C3BIC[3,1]

mu1.estimada.Grupo1.C3 = TabelaGrupo1C3BIC[1,2:dim(TabelaGrupo1C3BIC)[2]]
mu1.padronizada.Grupo1.C3 = mu1.estimada.Grupo1.C3/sum(TabelaGrupo1C3BIC[1,2:dim(TabelaGrupo1C3BIC)[2]])
mu2.estimada.Grupo1.C3 = TabelaGrupo1C3BIC[2,2:dim(TabelaGrupo1C3BIC)[2]]
mu2.padronizada.Grupo1.C3 = mu2.estimada.Grupo1.C3/sum(TabelaGrupo1C3BIC[2,2:dim(TabelaGrupo1C3BIC)[2]])
mu3.estimada.Grupo1.C3 = TabelaGrupo1C3BIC[3,2:dim(TabelaGrupo1C3BIC)[2]]
mu3.padronizada.Grupo1.C3 = mu3.estimada.Grupo1.C3/sum(TabelaGrupo1C3BIC[3,2:dim(TabelaGrupo1C3BIC)[2]])

for(i in 1:dim(U_ij[1:33,])[1]){
  somaGrupo1C3BIC[i,1] = 
    log(
      (theta1.Grupo1.C3 * sum((mu1.padronizada.Grupo1.C3^U_ij[i,]) ))+
        (theta2.Grupo1.C3 * sum((mu2.padronizada.Grupo1.C3^U_ij[i,]) ))+
        (theta3.Grupo1.C3 * sum((mu3.padronizada.Grupo1.C3^U_ij[i,]) ))
    )
}
somaGrupo1C3BIC
deltaC3Grupo1 = (2*kGrupo1)
BICGrupo1C3 = (2*sum(somaGrupo1C3BIC[,1]))-(deltaC3Grupo1*log(nGrupo1))
BICGrupo1C3

# Grupo 2
set.seed(13)
Grupo2C3inicialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[34:66,])[2]))
set.seed(1968)
Grupo2C3inicialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[34:66,])[2]))
set.seed(100)
Grupo2C3inicialC3 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[34:66,])[2]))

#Tempo grande para execução
Grupo2C3BIC = 
  SimulacoesC3EM(C = 3, nsim = 1, Matriz.real = U_ij[34:66,],
                 mu1.inicial = Grupo2C3inicialC1,
                 mu2.inicial = Grupo2C3inicialC2,
                 mu3.inicial = Grupo2C3inicialC3)
TabelaGrupo2C3BIC = cbind(rbind(mean(Grupo2C3BIC$theta[,1]),
                                mean(Grupo2C3BIC$theta[,2]),
                                mean(Grupo2C3BIC$theta[,3])),
                          Grupo2C3BIC$mu.media)
colnames(TabelaGrupo2C3BIC) = c("theta",dicionarioGrupo2);TabelaGrupo2C3BIC

somaGrupo2C3BIC = matrix(ncol = 1, nrow = dim(U_ij[34:66,])[1])

theta1.Grupo2.C3 = TabelaGrupo2C3BIC[1,1]
theta2.Grupo2.C3 = TabelaGrupo2C3BIC[2,1]
theta3.Grupo2.C3 = TabelaGrupo2C3BIC[3,1]

mu1.estimada.Grupo2.C3 = TabelaGrupo2C3BIC[1,2:dim(TabelaGrupo2C3BIC)[2]]
mu1.padronizada.Grupo2.C3 = mu1.estimada.Grupo2.C3/sum(TabelaGrupo2C3BIC[1,2:dim(TabelaGrupo2C3BIC)[2]])
mu2.estimada.Grupo2.C3 = TabelaGrupo2C3BIC[2,2:dim(TabelaGrupo2C3BIC)[2]]
mu2.padronizada.Grupo2.C3 = mu2.estimada.Grupo2.C3/sum(TabelaGrupo2C3BIC[2,2:dim(TabelaGrupo2C3BIC)[2]])
mu3.estimada.Grupo2.C3 = TabelaGrupo2C3BIC[3,2:dim(TabelaGrupo2C3BIC)[2]]
mu3.padronizada.Grupo2.C3 = mu3.estimada.Grupo2.C3/sum(TabelaGrupo2C3BIC[3,2:dim(TabelaGrupo2C3BIC)[2]])

for(i in 1:dim(U_ij[34:66,])[1]){
  somaGrupo2C3BIC[i,1] = 
    log(
      (theta1.Grupo2.C3 * sum((mu1.padronizada.Grupo2.C3^U_ij[i,]) ))+
        (theta2.Grupo2.C3 * sum((mu2.padronizada.Grupo2.C3^U_ij[i,]) ))+
        (theta3.Grupo2.C3 * sum((mu3.padronizada.Grupo2.C3^U_ij[i,]) ))
    )
}
somaGrupo2C3BIC
deltaC3Grupo2 = (2*kGrupo2)
BICGrupo2C3 = (2*sum(somaGrupo2C3BIC[,1]))-(deltaC3Grupo2*log(nGrupo2))
BICGrupo2C3

# Grupo 3
set.seed(13)
Grupo3C3inicialC1 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[67:100,])[2]))
set.seed(1968)
Grupo3C3inicialC2 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[67:100,])[2]))
set.seed(100)
Grupo3C3inicialC3 = rdirichlet(n = 1, alpha = rep(1,dim(U_ij[67:100,])[2]))

#Tempo grande para execução
Grupo3C3BIC = 
  SimulacoesC3EM(C = 3, nsim = 1, Matriz.real = U_ij[67:100,],
                 mu1.inicial = Grupo3C3inicialC1,
                 mu2.inicial = Grupo3C3inicialC2,
                 mu3.inicial = Grupo3C3inicialC3)
TabelaGrupo3C3BIC = cbind(rbind(mean(Grupo3C3BIC$theta[,1]),
                                mean(Grupo3C3BIC$theta[,2]),
                                mean(Grupo3C3BIC$theta[,3])),
                          Grupo3C3BIC$mu.media)
colnames(TabelaGrupo3C3BIC) = c("theta",dicionarioGrupo3);TabelaGrupo3C3BIC

somaGrupo3C3BIC = matrix(ncol = 1, nrow = dim(U_ij[67:100,])[1])

theta1.Grupo3.C3 = TabelaGrupo3C3BIC[1,1]
theta2.Grupo3.C3 = TabelaGrupo3C3BIC[2,1]
theta3.Grupo3.C3 = TabelaGrupo3C3BIC[3,1]

mu1.estimada.Grupo3.C3 = TabelaGrupo3C3BIC[1,2:dim(TabelaGrupo3C3BIC)[2]]
mu1.padronizada.Grupo3.C3 = mu1.estimada.Grupo3.C3/sum(TabelaGrupo3C3BIC[1,2:dim(TabelaGrupo3C3BIC)[2]])
mu2.estimada.Grupo3.C3 = TabelaGrupo3C3BIC[2,2:dim(TabelaGrupo3C3BIC)[2]]
mu2.padronizada.Grupo3.C3 = mu2.estimada.Grupo3.C3/sum(TabelaGrupo3C3BIC[2,2:dim(TabelaGrupo3C3BIC)[2]])
mu3.estimada.Grupo3.C3 = TabelaGrupo3C3BIC[3,2:dim(TabelaGrupo3C3BIC)[2]]
mu3.padronizada.Grupo3.C3 = mu3.estimada.Grupo3.C3/sum(TabelaGrupo3C3BIC[3,2:dim(TabelaGrupo3C3BIC)[2]])

for(i in 1:dim(U_ij[67:100,])[1]){
  somaGrupo3C3BIC[i,1] = 
    log(
      (theta1.Grupo3.C3 * sum((mu1.padronizada.Grupo3.C3^U_ij[i,]) ))+
        (theta2.Grupo3.C3 * sum((mu2.padronizada.Grupo3.C3^U_ij[i,]) ))+
        (theta3.Grupo3.C3 * sum((mu3.padronizada.Grupo3.C3^U_ij[i,]) ))
    )
}
somaGrupo3C3BIC
deltaC3Grupo3 = (2*kGrupo3)
BICGrupo3C3 = (2*sum(somaGrupo3C3BIC[,1]))-(deltaC3Grupo3*log(nGrupo3))
BICGrupo3C3

####============
#### Resultados
####============
rbind(BICGrupo1C1,BICGrupo1C2,BICGrupo1C3)
rbind(BICGrupo2C1,BICGrupo2C2,BICGrupo2C3)
rbind(BICGrupo3C1,BICGrupo3C2,BICGrupo3C3)

resultados_aloc_EM = AlocLikelihoodIndC3(MatrizU = U_ij,
                                         Thetas = c(0.33,0.33,0.34),
                                         mu1.estimate = medida1,
                                         mu2.estimate = medida2,
                                         mu3.estimate = medida3)
