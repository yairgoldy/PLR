### R/PLR: Estimation in the Proportional Likelihood Model with Missing Data and Biased Sampling.



[Yair Goldberg](https://yairgo.net.technion.ac.il/) and [Malka Gorfine](https://yairgo.net.technion.ac.il/)

[R/PLR](https://github.com/yairgoldy/PLR) is an [R](https:/www.r-project.org) package. This package implements estimation of the parameters of the proportional likelihood model with Lou & Tsai (2012), Chan (2013), and Goldberg & Gorfine (2019) methods. It can handle missing data and biased sampling. The package also contains a function that creates the different datasets used in the simulation of oldberg & Gorfine (2019).


#### Installation

You can install it from its [GitHub repository](https://github.com/yairgoldy/PLR). You first need to install the [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install R/corihw using the `install_github` function in the
[devtools](https://github.com/hadley/devtools) package.

```r
library(devtools)
install_github("yairgoldy/PLR")
```

#### Example use

Create data from Setting 4 of Goldberg & Gorfine (2019) and calculate the estimator of $\beta$ using  Lou & Tsai (2012), Chan (2013), and Goldberg & Gorfine (2019) methods. Here, the covariate vector consists of $X=(X_1,X_2)^T$, where $X_2$ follows a zero-mean normal distribution with standard deviation 0.5, and given $X_2$, $X_1$ follows Bernoulli distribution with success probability $\exp(1-X_2)/\{1+\exp(1-X_2)\}$. The true parameters' values are $\beta=(\beta_1,\beta_2)^T=(-1,-1)^T$. $Y$ is continuous and the baseline density is defined by 
$$
g(y) = \{\Phi(0.5)\}^{-1} (0.5 \pi)^{-1/2} \exp\{-2(y-0.25)^2\} \;\;\;\;\;\; y \geq 0  \, , 
$$
where $\Phi$ is the standard normal cumulative distribution function. The probability of observing complete data consists of $$\mbox{pr}(R=1|Y,X_1,X_2)=\frac{\exp(1-X_2-2Y)}{1+\exp(1-X_2-2Y)}\,,$$

The following code creates 50 dataset, each of sample-size 200 and compare the performence of the three methods.



```{r}

library(PLR)
library(tidyverse)


n <- 200
nsim <- 50

# Create an empty table to hold the results
res <- tibble(method = character(),  Beta1 = numeric(), Beta2 = numeric())


for(seed in 1:nsim){
  set.seed(seed)
  dat <- generate_example(n = n, setting = 4,
                          parameters = 0.5, seed = seed)
  for(method in c("Lou.Tsai","Chan","Gold.Gorf")){
    
    # Compute the estimate of beta and the baseline distribution G
    est <-  lpr(x = dat$x,y = dat$y, method = method, weights = dat$weights)
    beta <- est$beta
    
    # Add the results to the table
    res <- rbind(res, tibble(method = method,  Beta1 = beta[1], Beta2 = beta[2]))
  }
}
  res <- gather(res,Parameter,Estimate,-method)


```

Draw a graph of the results
```{r}

ggplot(res,aes(y=Estimate,fill=method))+geom_boxplot()+facet_wrap(~Parameter)+geom_hline(yintercept = -1)

```




