plot_pcorgraph <- function(Theta, pos_clr, neg_clr, plot_layout, label_cex){

  K <- length(Theta)
  
  A_list <- Theta
  G_list <- vector("list", K)
  
  for(k in 1:K){
    A_list[[k]][which(A_list[[k]] != 0)] <- 1
    G_list[[k]] <- graph_from_adjacency_matrix(A_list[[k]], mode = "undirected", diag = FALSE)
  }
  
  pcor_full <- vector("list", K)
  for(k in 1:K){
    pcor_full[[k]] <- Theta[[k]]
  }
  
  # obtain partial correlations
  scaled_pcor_full <- pcor_full
  for(k in 1:K){
    for(i in 1:nrow(pcor_full[[k]])){
      for(j in 1:ncol(pcor_full[[k]])){
        scaled_pcor_full[[k]][i,j] <- -1*(pcor_full[[k]][i,j]/(sqrt(diag(pcor_full[[k]])[i])*sqrt(diag(pcor_full[[k]])[j])))
      }
    }
    diag(scaled_pcor_full[[k]]) <- diag(pcor_full[[k]])
  }
  
  # edge width weights
  for(k in 1:K){
    E(G_list[[k]])$e.weight <- 10*abs(scaled_pcor_full[[k]][lower.tri(scaled_pcor_full[[k]])][which(scaled_pcor_full[[k]][lower.tri(scaled_pcor_full[[k]])] !=0)])
  }
  
  # edge color weights
  for(k in 1:K){
    E(G_list[[k]])$ec.weight <- scaled_pcor_full[[k]][lower.tri(scaled_pcor_full[[k]])][which(scaled_pcor_full[[k]][lower.tri(scaled_pcor_full[[k]])] !=0)]
  }
  
  # assign color
  for(k in 1:K){
    for(i in 1:length(E(G_list[[k]])$ec.weight)){
      if(E(G_list[[k]])$ec.weight[i] >= 0) E(G_list[[k]])$color[i] <- pos_clr else E(G_list[[k]])$color[i] <- neg_clr
    }
  }
  
  # plot graph
  coords <- layout.fruchterman.reingold(G_list[[1]])
  par(mfrow = plot_layout)
  
  for(k in 1:K){
    plot(G_list[[k]],vertex.color="white", vertex.size=degree(G_list[[k]])*.35, vertex.label.cex = label_cex, vertex.label.font = 1, layout = coords, edge.curved = 0.5, edge.color = E(G_list[[k]])$color, edge.width = E(G_list[[k]])$e.weight, vertex.label.color= "black", vertex.label.dist=0.5, vertex.label.degree=pi)
  }
}
