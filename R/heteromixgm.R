heteromixgm <- function(X, method, lambda1, lambda2, ncores){
  if(method == "Gibbs")
  {
    K <- length(X)
    if(length(lambda1) == 1 & length(lambda2) == 1)
    {
      K <- length(X)
      Theta <- vector(mode = "list", length = K)
      for(k in 1:K)
      {
        Theta[[k]] <- diag(ncol(X[[k]]))
        Theta[[k]] <- as(Theta[[k]], "sparseMatrix")
      }
      est = Gibbs_method(X, lambdas=c(lambda1, lambda2), n_lambda = 1, Theta = Theta, combination = 1, max.elongation = 5, em.tol=0.001, ncores = ncores)
    }
    else
    {
      lambdas <- expand.grid(lambda1, lambda2)
      n_lambda <- nrow(lambdas)
      est = vector("list", n_lambda)
      for(combination in 1:n_lambda)
      {
        if(combination == 1)
        {
          est[[combination]] = vector("list", n_lambda)
          Theta <- vector(mode = "list", length = K)
          for(k in 1:K){
            Theta[[k]] <- diag(ncol(X[[k]]))
            Theta[[k]] <- as(Theta[[k]], "sparseMatrix")
          }
          m <- paste(c("Obtaining parameter estimates using the Gibbs method:", floor(100 * combination/n_lambda), "%", "\n"), collapse="")
          cat(m, "\r")
          est[[combination]] = Gibbs_method(X, lambdas=lambdas, n_lambda=n_lambda, Theta = Theta, combination = combination, max.elongation = 5, em.tol=0.001, ncores = ncores)
        }else{
          est[[combination]] = vector("list", n_lambda)
          for(k in 1:K){
            Theta[[k]] = est[[(combination - 1)]]$Theta[[k]]
            Theta[[k]] = as(Theta[[k]], "dgTMatrix")
            Theta[[k]] = as(Theta[[k]], "sparseMatrix")
          }
          m <- paste(c("Obtaining parameter estimates using the Gibbs method:", floor(100 * combination/n_lambda), "%", "\n"), collapse="")
          cat(m, "\r")
          est[[combination]] = Gibbs_method(X, lambdas=lambdas, n_lambda=n_lambda, Theta= Theta, combination = combination, max.elongation = 5, em.tol=0.001, ncores = ncores)
        }
      }
    }
  }
  else if(method == "Approximate")
  {
    if(length(lambda1) == 1 & length(lambda2) == 1)
    {
      ini <- initialize(X, ncores = ncores)
      Z	= ini$Z
      ES	= ini$ES
      lower_upper = ini$lower_upper
      est = Approx_method(X, Z, ES=ES, lambdas=c(lambda1, lambda2), lower_upper=lower_upper, combination = 1, em_tol=0.001, em_iter=5, ncores = ncores)
    }
    else
    {
      lambdas <- expand.grid(lambda1, lambda2)
      n_lambda <- nrow(lambdas)
      ini <- initialize(X, ncores = ncores)
      Z	= ini$Z
      ES	= ini$ES
      lower_upper = ini$lower_upper

      est = vector("list", n_lambda)
      for(combination in 1:n_lambda)
      {
        m <- paste(c("Obtaining parameter estimates using the approximate method:", floor(100 * combination/n_lambda), "%", "\n"), collapse="")
        cat(m, "\r")
        est[[combination]] = vector("list", n_lambda)
        est[[combination]] = Approx_method(X, Z, ES=ES, lambdas=lambdas, lower_upper=lower_upper, combination = combination, em_tol=0.001, em_iter=5, ncores = ncores)
      }
    }
  }
  else
  {
    message("Please use either Gibbs or Approximate method")
  }
  return(est)
}
