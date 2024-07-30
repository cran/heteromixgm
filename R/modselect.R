modselect <- function(est, X, l1, l2, gamma){

  K <- length(X)
  n_lambda <- nrow(expand.grid(l1, l2))
  thetas <- vector("list", n_lambda)
  A <- vector("list", n_lambda)

  K <- length(X)
  p <- ncol(X[[1]])
  n <- sapply(X, nrow)

  nu <- matrix(NA, n_lambda, K)
  lik_mat <- matrix(NA, n_lambda, K)
  lik <- numeric(n_lambda)
  aic <- numeric(n_lambda)
  bic <- numeric(n_lambda)
  ebic <- numeric(n_lambda)

  lognu <- matrix(nrow = n_lambda, ncol = K)
  rhs <- matrix(nrow = n_lambda, ncol = K)

  for(i in 1:n_lambda)
  {
    for(k in 1:K){
      thetas[[i]][[k]] <- as.matrix(est[[i]]$Theta[[k]])
      A[[i]][[k]] <- thetas[[i]][[k]]
      # adjacency matrix
      A[[i]][[k]][which(A[[i]][[k]] != 0)] <- 1
      # nonzero elements penalty
      nu[i,k] <- sum(A[[i]][[k]][upper.tri(A[[i]][[k]])])
      # likelihood
      lik_mat[i,k] <- -n[k] *(determinant(as.matrix(thetas[[i]][[k]]), logarithm = TRUE)$modulus - sum(diag(as.matrix(est[[1]]$ES[[k]]) %*%as.matrix(thetas[[i]][[k]]))))

      # (e)BIC penalties
      lognu[i,k] <- log(n[k]) * nu[i,k]
      rhs[i,k] <- (4 * gamma * log(p) * nu[i,k])
    }

    lik[i] <- sum(lik_mat[i,])

    # information criteria
    aic[i] <- ( lik[i] ) + ( 2 * sum(nu[i,]) )
    bic[i] <- lik[i] + sum(lognu[i,])
    ebic[i] <- lik[i] + sum(lognu[i,]) + sum(rhs[i,])
  }

  theta_aic <- est[[which.min(aic)]]$Theta 
  theta_bic <- est[[which.min(bic)]]$Theta 
  theta_ebic <- est[[which.min(ebic)]]$Theta 

  returnlist <- list("theta_aic" = theta_aic, "theta_bic" = theta_bic, "theta_ebic" = theta_ebic)
  return(returnlist)
}

