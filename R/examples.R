

#' 'generate_example': Generate the data examples presented in Goldberg and Gorfine (2019)
#' @param n  Sample size
#' @param beta  A 2x1 vector of coefficeients.
#' @param setting  The setting number is a number between 1 and 4.
#' See the Simulation Section of Goldberg and Gorfine (2019) for the description of the settings.
#' @param parameters The parameters of the example. The parameters depend on the example's setting.
#' For setting 1,3, and 4 the parameter that was used is 0.25.
#' For setting 2, the parameter that was used is 3
#' @param seed Use the seed in creating the example. If seed is negetive, no seed is used.
#' @return A list that contains the nx2 matrix of explanatory variables x, the response vector y.
#' If applicable also the missingness indicator R and a vector of weights.
#' @importFrom tibble tibble
#' @importFrom stats rbinom rexp
#' @export
generate_example<- function(n = 100, beta= c(-1,-1), setting = 1, parameters = 1, seed =-1){
  if(seed>0){  set.seed(seed = seed)}
  z <- generate_z(n)
  y <- generate_y(n,beta,parameters[1],z,setting)
  if(setting == 1 || setting == 2){
    return(list(x = z,y=y))
  } else if(setting == 3){
    z2 <- z[,2]
    eta <- 1-z2
    R <- rbinom(n = n,size=1,exp(eta)/(1+exp(eta)))


    dat <- cbind(z,y,R)
    R <- dat[,4]
    x <- dat[R==1,1:2]
    z2 <- dat[,2]
    y <- dat[R==1,3]

    #     dat <- tibble(z2,R)

    logistic_model <- glm(R~z2,family = binomial(link = "logit"))
    ps <- predict(logistic_model,type="response",newdata = tibble(z2))

    return(list(x = x,y = y, R = R[R==1],weights=1/(ps[R==1])))

  } else if(setting == 4){
    z2 <- z[,2]
    eta <- 1-z2-2*y

    R <- rbinom(n = n,size=1,exp(eta)/(1+exp(eta)))


    dat <- cbind(z,y,R)
    R <- dat[,4]
    x <- dat[R==1,1:2]
    z2 <- dat[,2]
    y <- dat[,3]

    #     dat <- tibble(z2,R)

    logistic_model <- glm(R~z2+y,family = binomial(link = "logit"))
    ps <- predict(logistic_model,type="response",newdata = tibble(z2,y))
    y <- dat[R==1,3]
    return(list(x = x,y = y, R = R[R==1],weights=1/(ps[R==1])))

  }
}




generate_z <- function(n){
  z2 <- rnorm(n,mean=0,sd=0.5);
  prob <- exp(1-z2)/(1+exp(1-z2));
  z1 <- rbinom(n,size=1,prob=prob);
  cntr_z <- cbind(z1,z2);
  return(cntr_z);
}

generate_y <- function(n,beta,parameter,z,setting){
  if(setting == 2){
    y <- numeric(n)
    for(i in 1:n){
      probs.y<- gy(beta,z[i,],parameter)
      y[i] <- sample(0:100,size=1,prob=probs.y)
    }
  } else {
    u <- runif(n)
    c <- z%*%beta;
    y <- numeric(n)
    for(i in 1:n){
      sol<- uniroot(eval_y,lower=0,upper=30,c=c[i],s2=parameter,u=u[i]);
      y[i] <- sol$root
    }
  }
  return(y)
}



erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}

CDF_Y_given_Z <- function(c,s2,y){
  -(erf((c*s2+s2-y)/sqrt(2*s2))-erf((c*s2+s2-0)/sqrt(2*s2)))/(1+erf((c*s2+s2-0)/sqrt(2*s2)))
}
eval_y <- function(y,c,s2,u){
  CDF_Y_given_Z(c,s2,y)-u;
}



gy<- function(beta,z,lambda){
  y <- 0:100
  c <- sum(z*beta);
  return((1+y)*(lambda^y)*exp(c*y)/
           ( factorial(y)*exp(lambda*exp(c)) * (1+lambda*exp(c) )) )
}


gdist <- function(gfun,setting){
  if(setting == 2){
    s <- 0:100
    sum( abs (gfun(s)-cumsum(gy(c(-1,-1),c(0,0),3))  )  )

  } else {
    f <- function(y){abs (gfun(y)-CDF_Y_given_Z(0,0.25,y))}
    return(integrate(f, 0, 5,rel.tol=.Machine$double.eps^.05))
  }
}
