####============================================
#### Aplicação da modularidade - César Macieira
####============================================
rm(list=ls(all=T))

# Packages
if(!require(igraph)){ install.packages("igraph"); require(igraph)}

# 1. Inicialização
set.seed(123)
n_pessoas <- 10
n_palavras <- 9
matriz_medidas <- matrix(0, nrow = n_pessoas, ncol = n_palavras)

# 2. Gerar as medidas por grupo
for (i in 1:n_pessoas) {
  if (i <= 3) {
    matriz_medidas[i, 1:3] <- runif(3)
  } else if (i <= 6) {
    matriz_medidas[i, 4:6] <- runif(3)
  } else {
    matriz_medidas[i, 7:9] <- runif(3)
  }
}

# 3. Criação da matriz U
matriz_binaria <- matrix(0, nrow = n_pessoas, ncol = n_palavras)
for (i in 1:n_pessoas) {
  top_3 <- order(matriz_medidas[i, ], decreasing = TRUE)[1:3]
  matriz_binaria[i, top_3] <- 1
}
matriz_binaria

# 4. Criar grafo bipartido
linhas <- which(matriz_binaria == 1, arr.ind = TRUE)
edges <- cbind(paste0("p", linhas[,1]), paste0("w", linhas[,2]))
g_bipartido <- graph_from_edgelist(edges, directed = FALSE)
V(g_bipartido)$type <- grepl("^w", V(g_bipartido)$name)

# 5. Projetar grafo entre pessoas
g_pessoas <- bipartite_projection(g_bipartido)$proj1

# 6. Comunidades e modularidade
comunidades <- cluster_louvain(g_pessoas)
modularidade <- modularity(comunidades)
cat("Modularidade:", modularidade, "\n")

# 7. Plot do grafo
plot(
  comunidades, 
  g_pessoas,
  vertex.label = V(g_pessoas)$name,
  vertex.color = comunidades$membership,
  vertex.size = 30,
  edge.width = 2,
  main = "Grafo de Pessoas com Comunidades (Louvain)"
)
