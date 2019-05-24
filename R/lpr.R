
#' lpr: Main function for calculating the model parameters for a proportional likelihood model
#'
#'
#'
#' @param x  An nxp matrix of the explanatory variables
#' @param y  An nx1 vector of response
#' @param method The method to use when calculating the parameters:
#'  "1" or "Lou.Tsai" for the method of Lou & Tsai (2012),
#'  "2" or "Chan" for the method of Chen (2013),
#'  "3" or "Gold.Gorf" for the method of Goldberg and Gorfine (2019).
#'  Method of Gold.Gorf can handle survival. Methods Chan and Gold.Gorf can handle missing data.
#' @param weights Vector of weights for missing data or survival. Can be used only with Gold.Gorf method.
#' The weights should be the inverse of the probability of missingness or censoring.
#'
#' @return A list that contains a px1 vector beta and the baseline distribution G (see Lou & Tsai, 2012, for details)
#'
#' @importFrom dplyr group_by select filter mutate summarise arrange ungroup
#' @importFrom tibble tibble
#' @import tidyr
#' @importFrom stats approxfun optim binomial glm integrate optimise pnorm predict rnorm runif uniroot
#' @export
lpr <- function(x,y,method,weights=NULL){

  p <- ncol(x)
  n <- nrow(x)

  if(is.null(weights)) weights = rep(1,n)


  if( ! is.matrix(x)){
    x <- matrix(x,ncol = 1)
  }
  switch (as.character(method),
          "1" =,
          "Lou.Tsai" = {fun <-  lt_fun;freq.list <- freq(x,y);weights = rep(1,n)},
          "2" =,
          "Chan" = {fun <-  chan_fun; freq.list <- freq(x,y);weights = rep(1,n)},
          "3" =,
          "Gold.Gorf" = {
            freq.list <- freq(x,y)
            if(length(freq.list$yk)<50){
              fun <-  gg_discrete_fun
            } else {
              fun <-  gg_continuous_fun
            }
          },
          stop("Please choose a method from Lou.Tsai, Chan, and Gold.Gorf")
  )



  if(p == 1){
    sol <- optimise(f = fun,interval = c(-10,10),x=x,y=y,freq.list=freq.list,weights = weights,list(maxit = 50))
    beta <- sol$minimum
  } else {
    sol <- optim(rep(0,p),fun,x=x,y=y,freq.list=freq.list,weights = weights)
    beta <- sol$par
  }

  g <- profile.p.gg(y,x,beta,weights)
  g_fun <- tibble(y,g) %>%
    group_by(y) %>%
    summarise(g=sum(g)) %>%
    arrange(y) %>%
    mutate(g=cumsum(g))


  g <- approxfun(g_fun$y, g_fun$g,method="constant", rule = 2)
  return(list(beta = beta, g = g))

}
