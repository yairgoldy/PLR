lt_fun <- function(beta,x,y,freq.list,weights){
  zl <- freq.list$zl
  lambdal <- freq.list$lambdal
  yk <- freq.list$yk
  nk <- freq.list$nk
  nkl <- freq.list$nkl
  p <- profile.p(nk,yk,zl,lambdal,beta)
  
  p.loglike <- 0
  for(i in 1:nrow(zl)){
    p.loglike <- p.loglike+sum(nkl[,i]*yk)*sum(zl[i,]*beta)-lambdal[i]*log(sum(p*exp(sum(zl[i,]*beta)*yk)))
  }
  p.loglike <- p.loglike+sum(nk*log(p))
  return(-p.loglike)
}


chan_fun <- function(beta,x,y,freq.list,weights){
  
  res <- numeric(length(beta))
  for(i in 1:length(y)){
    for(j in 1:length(y)){
      s <- (x[i,]-x[j,])*(y[i]-y[j])
      res <- res+s* exp(-sum(s*beta))/(1+exp(-sum(beta*s)))
    }
  }
  
  
  return( sum(res^2))
}


gg_discrete_fun <- function(beta,x,y,freq.list,weights)
{ 
  n <- length(y)
 
  p <- ncol(x)
  S0 <- numeric(n)
  EYgivenX <- numeric(n)
  alpha <- matrix(0,n,p)
  
  zl <- freq.list$zl
  lambdal <- freq.list$lambdal
  yk <- freq.list$yk
  nk <- freq.list$nk
  nkl <- freq.list$nkl
  pk <- profile.p(nk,yk,zl,lambdal,beta)
  py <- numeric(n)
  K <- length(yk)
  
  
  U0 <- numeric(length(yk))
  U1 <- matrix(0,length(yk),p)
  V1 <- matrix(0,n,p)
  weights_k <- numeric(K)
  
  
  for(j in 1:K){
    y.indeces<-abs(y-yk[j])<0.0001 
    py[y.indeces] <- pk[j]/nk[j]   
    weights_k[j] <- sum(weights[y.indeces == T])
  }
  
  
  
  
  
  for(j in 1:n){
    y.indeces<-abs(y-yk[j])<0.0001 
    py[y.indeces] <- pk[j]/nk[j]    
  }
  
  
  disc.mat <- matrix(0,n,K)
  for (i in 1:n) {
    for (j in 1:K) {
      disc.mat[i,j] <- exp(sum(beta*x[i,]*yk[j]))*pk[j] 
    }
  }
  for(i in 1:n){
    
    S0[i] <- sum(  py*exp( sum(beta*x[i,])*y   ))
    EYgivenX[i] <- sum(py*y*  exp( sum(beta*x[i,])*y) )/S0[i]
    alpha[i,] <- x[i,]*(y[i]- EYgivenX[i]  )
  }
  for(j in 1:K){ 
    disc <- numeric(n)
    for(i in 1:n) {
      disc[i] <- 1 - disc.mat[i,j] / sum(disc.mat[i,j:K]) 
      if (j > 1) disc[i] <- disc[i] + disc.mat[i,j-1] / sum(disc.mat[i,(j-1):K])
    }
    U0[j] <- mean(weights*disc* as.vector( exp( (x%*%beta) * yk[j] )/S0))
    U1[j,] <- apply(weights*  x* ( ( rep(yk[j],n)-EYgivenX)*disc*as.vector(exp(  (x%*%beta)* yk[j]))/S0),2, mean)
  }
  for(i in 1:n){
    V1[i,] <- (1/S0[i]) *apply( (U1/U0)* as.vector(pk*exp( sum(beta*x[i,])*yk)),2,sum)
  }
  
  return( sum((apply( weights*(alpha+V1),2,sum) - apply(weights_k*U1/U0,2,sum))^2));
}



gg_continuous_fun <-  function(beta,x,y,freq.list,weights){ 
  n <- length(y)
  p <- ncol(x)
  S0 <- numeric(n)
  EYgivenX <- numeric(n)
  alpha <- matrix(0,n,p)
  U0 <- numeric(n)
  U1 <- matrix(0,n,p)
  V1 <- matrix(0,n,p)
  
  py <- profile.p.gg(y,x,beta,weights)
  
 
  for(i in 1:n){
    
    S0[i] <- sum( weights*py*exp( sum(beta*x[i,])*y   ))
    EYgivenX[i] <- sum(weights*py*y*  exp( sum(beta*x[i,])*y) )/S0[i]
    alpha[i,] <- x[i,]*(y[i]- EYgivenX[i]  )
  }
  for(i in 1:n){
    U0[i] <- mean( weights*exp( (x%*%beta) * y[i] )/S0)
    U1[i,] <- apply( weights* x* ( ( rep(y[i],n)-EYgivenX)*as.vector(exp(  (x%*%beta)* y[i]))/S0),2, mean)
  }
  for(i in 1:n){
    V1[i,] <- (1/S0[i]) *apply( weights*(U1/U0)* as.vector(py*exp( sum(beta*x[i,])*y)),2,sum)
  }
  return( sum((apply( (alpha- U1/U0+V1)*weights,2,mean))^2));
}



