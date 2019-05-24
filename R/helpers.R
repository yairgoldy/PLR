
freq <- function(z,y){
  ## Compute zl and lambdal

  # 1. We order the matrix by a lexicographic order
  ordered.zy  <- ord.mat(cbind(z,y),cols = 1:ncol(z) )
  orderedz <- matrix(ordered.zy[,1:ncol(z)],ncol=ncol(z))
  orderedy <- ordered.zy[,ncol(z)+1]


  # 2. We run through the rows of the matrix
  #    a. we put all unique values in zl
  #    b. we put the number of repeatition of the unqiue values in lambdal
  #    c. we create the vector names_zl of length n which assigns the
  #       number to each unique row
  previousz <- orderedz[1,]
  lambdal <- 1
  zl <- orderedz[1,]
  names.zl <-numeric(nrow(z))
  names.zl[1] <- 1

  for(i in 2:nrow(orderedz)){
    currentz <- orderedz[i,]
    if( norm(as.matrix(currentz-previousz))==0){
      lambdal[length(lambdal)] <- lambdal[length(lambdal)]+1
      names.zl[i] <-  names.zl[i-1]
    }
    else{
      zl <- rbind(zl,currentz);
      lambdal <- c(lambdal,1);
      names.zl[i] <-  names.zl[i-1]+1
    }
    previousz <- currentz;
  }



  ## Compute yk and nk
  ty <- table(orderedy)
  nk <- as.vector(ty)
  yk <-as.numeric(names(ty))


  nkl <- as.matrix(table(orderedy,names.zl))
  if( !is.matrix(zl)){zl <- matrix(zl,ncol = 1)}
  return(list(zl=zl,lambdal=lambdal,yk=yk,nk=nk,nkl=nkl))
}


Afun <- function(y,x,beta,p,i,weights=1){
  apply(x,1,function(xi){sum(weights*y^i * p*exp(sum(beta*xi)*y)) })
}

profile.p <- function(nk,yk,zl,lambdal,beta){
  condition <- T
  p.old <- nk/sum(nk)
  num.itr <- 0
  while(condition==T){

    a0 <- Afun(yk,zl,beta,p.old,0)

    p.new <- nk/ apply(as.matrix(yk),1,function(yi){ sum(lambdal*exp( (zl%*%as.matrix(beta,ncol=1))*yi)/a0 )})
    p.new <- p.new/sum(p.new)
    num.itr <- num.itr+1
    condition <- (max(abs(p.new-p.old))>0.00001) & num.itr<100
    p.old <- p.new
  }
  return(p.new)
}

profile.p.gg <- function(y,x,beta,weights=1){
  condition <- T
  n <- length(y)

  p.old <- rep(1/n,n)
  num.itr <- 0
  while(condition==T){

    a0 <- Afun(y,x,beta,p.old,0,weights)

    p.new <- weights/ apply(matrix(y,ncol = 1),1,function(yi){ sum (weights*exp( (x%*%matrix(beta,ncol=1))*yi)/a0 )})
    p.new <- p.new/sum(p.new)
    num.itr <- num.itr+1
    condition <- (max(abs(p.new-p.old))>0.00001) & num.itr<100
    p.old <- p.new
  }

  return(p.new)
}




# Afun <- function(yk,zl,beta,p,i){
#   apply(zl,1,function(zi){sum(yk^i * p*exp(sum(beta*zi)*yk)) })
# }
#
#
# profile.p <- function(nk,yk,zl,lambdal,beta){
#   condition <- T
#   p.old <- nk/sum(nk)
#   num.itr <- 0
#   while(condition==T){
#
#     a0 <- Afun(yk,zl,beta,p.old,0)
#
#     p.new <- nk/ apply(as.matrix(yk),1,function(yi){ sum(lambdal*exp( (zl%*%as.matrix(beta,ncol=1))*yi)/a0 )})
#     p.new <- p.new/sum(p.new)
#     num.itr <- num.itr+1
#     condition <- (max(abs(p.new-p.old))>0.00001) & num.itr<100
#     p.old <- p.new
#   }
#   return(p.new)
# }


ord.mat  <-  function(M, decr = F, cols = NULL){
  if(is.null(cols))
    cols = 1: ncol(M)
  out = do.call( "order", as.data.frame(M[,cols]))
  if (decr)
    out = rev(out)
  return(M[out,])
}

freq.censored <- function(z,y,delta){
  ## Compute zl and lambdal

  # 1. We order the matrix by a lexicographic order
  m <- cbind(z,y)
  ordered.zy  <- m[order(y), ]
  zl <- ordered.zy[,1:ncol(z)]
  lambdal <- rep(1,nrow(z))
  orderedy <- ordered.zy[,ncol(z)+1]





  ## Compute yk and nk
  ty <- table(orderedy)
  nk <- as.vector(ty)
  yk <-as.numeric(names(ty))


  nkl <- matrix(0,nrow=length(yk),ncol=nrow(z))
  mkl <- matrix(0,nrow=length(yk),ncol=nrow(z))
  for(k in 1:length(yk)){
    nkl[k,]  <- (abs(y-yk[k])< 0.00001 & delta==T)
    mkl[k,]  <- (abs(y-yk[k])< 0.00001 & delta==F)

  }
  if( !is.matrix(zl)){zl <- matrix(zl,ncol = 1)}

  return(list(zl=zl,lambdal=lambdal,yk=yk,nk=nk,nkl=nkl,mkl=mkl))
}



profile.p.zhu <- function(nk,yk,zl,lambdal,beta,nkl,mkl){
  condition <- T
  p.old <- nk/sum(nk)
  num.itr <- 0
  while(condition==T){

    a0 <- Afun(yk,zl,beta,p.old,0)
    # this are the first and second term in the denominater on page page 2473, Zhu 2014
    expr1 <- apply(as.matrix(yk),1,function(yi){ sum(lambdal*exp( (zl%*%as.matrix(beta,ncol=1))*yi)/a0 )})
    v <- numeric(length(yk))



    for(i in 1:nrow(zl)){
      denumer <-  rev(cumsum(rev(p.old* exp(sum(zl[i,]*beta) *yk) )))
      numer <- mkl[,i]*exp(sum(zl[i,]*beta)*yk)
      v <- v+numer/denumer
    }
    expr2 <- sum(v)
    p.new <- nk/(expr1-expr2)
    p.new <- max(p.new,0.00000001)
    p.new <- p.new/sum(p.new)
    num.itr <- num.itr+1
    condition <- (max(abs(p.new-p.old))>0.00001) & num.itr<100
    p.old <- p.new

  }

  return(p.new)
}




