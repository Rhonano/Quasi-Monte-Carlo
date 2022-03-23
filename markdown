---
title: 
authors:
  - name: 
    thanks: 
    department: 
    affiliation: 
    location: 
    email: 
  - name: 
    department: 
    affiliation: 
    location: 
    email: 
abstract: |
  Enter the text of your abstract here.
bibliography: references.bib
biblio-style: unsrt
output: 
  pdf_document: 
    number_sections: true
toc: true
---
```{r include=FALSE}
library(fOptions)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(LSMonteCarlo)
library(Jdmbs)
library(knitr)
library(knitr)
library(bookdown)
library(gridExtra)
```

# Introduction

Here goes an introduction text

# American option valuation problem
In order to understand how to valuate American put options it is necessary to define how it differs from the more simple European options. A put option is a security that gives the owner the right to sell a fixed number of shares of an underlying stock at a fixed price at any time before or on the given maturity date. There are many types of options where the most common two are European and American. The European type can only be exercised on the maturity date whereas the American type can be exercised at any time up until the expiration. The pricing problem of American options can be simplified to specifying a process $U(t), 0\leq t \leq T$ which represents the discounted payoff from exercise at time t, and a class of stopping times $\mathcal T$. The problem is then to find the optimal expected discounted payoff$$\sup_{\tau\in \mathcal T}E[U(\tau)]$$

# Least squares Monte Carlo

## Least Square Monte Carlo algorithm
The LSM approach uses least squares to recursively approximate the discounted payoff at $t_{k-1},t_{k-2},...,t_1$. The reason for working backwards is that the cash flows can differ since it may be optimal to stop before the maturity date. We assume that the unknown functional form for approximating conditional expectation can be represented as a linear combination of basis functions. Assuming that the value of the asset underlying the option follows a Markov process^[A Markov chain is a stochastic model describing a sequence of possible events in which the probability of each event depends only on the state attained in the previous event (Gagniuc, 2017)]. 

## A numerical example
The key idea of the least squares Monte Carlo (LSM) approach is to estimate the conditional expected discounted payoff from the cross-sectional information of the Monte Carlo simulation. 

For a better understanding I produce a simple example below. Consider an American put option on a share of a non-dividend-paying stock. The option is exercisable at times t=1,2,3 where 3 is the maturity date. The paths are generated using a geometric Brownian motion with a risk free rate of 6%. 

```{r fig1, echo=FALSE}
set.seed(2)
n <- 8
GBM <- rnorm.pseudo(n = 8, dimension = 3)
for (i in 1:n) {
  GBM[i, ] <- 1 * exp(cumsum(((0.1 - 0) * (1/365) - 0.5 * 
                                   0.4 * 0.4 * (1/365)) + (0.4 * (sqrt(1/365)) * 
                                                                GBM[i, ])))}

sim <- as.data.frame(GBM)
sim <- rename(sim, "t=1"="V1","t=2"="V2","t=3"="V3")
num_ex <- round(cbind(as.data.frame(matrix(1,8,1)),sim),3)
num_ex <- rename(num_ex, "t=0"="V1")

kable(num_ex, caption="Stock price paths")
```
The objective is then to recursively compute intermediate matrices to solve for the stopping rule that maximizes the value of the option at each time t, beginning at t=3. These cash flows are also equal to that of the European option. 


```{r echo=FALSE}
empty <- rename(as.data.frame(matrix("--",8,2)),"t=1"="V1","t=2"="V2")
cflow <- pmax(0,1.05-num_ex["t=3"])
cflow3 <- cbind(empty, cflow)
cflow3[is.na(cflow3)]=0
cflow3 <- rename(cflow3,"t=3"="cflow")
kable(cflow3, caption="Cash flow at t=3")
```
The 0 values indicate that the path is out-of-money and not worth anything. By only using the paths that are in-the-money I get a better estimate of the conditional expectation function in addition to speeding up the algorithm. Let X denote the stock prices as time 2 and Y denote the discounted cash flows at time 3 if not exercised at time 2. 

```{r echo=FALSE}
Y <- round(cflow3["t=3"]/(1+(0.1)/3),3)
X <- num_ex["t=2"]
X[3,1] <- "--"
X[7,1] <- "--"
reg2 <- cbind(Y,X)
reg2 <- rename(reg2,"Y"="t=3","X"="t=2")

kable(reg2,caption="Regression at time 2")

```

I can then regress Y on a constant X and X^2 which in this case results in the function $E[Y|X]=-13.35+27.98X-14.59X^2$. With this simple conditional expectation model it is possible to compare the value of immediate exercise at time 2, with the value from continuation. 
```{r echo=FALSE}
fit_data <- reg2[-c(3,7),]
fit_data <- data.frame(lapply(fit_data,as.numeric))
fit_data <- fit_data %>% 
  mutate("XX"=X^2)
q_model <- lm(formula = Y~X+XX, data=fit_data)

ex_2 <- round(pmax(0,1.05-num_ex["t=2"]),3)
ex_2[3]="--"
ex_2[7]="--"
pred <- c(0.006,0.040,"--",0.044,0.051,0.058,"--",0.063)
opt_early <- as.data.frame(cbind(ex_2,pred))
opt_early <- rename(opt_early,"Exercise"="ex_2","Continuation"="pred")
kable(opt_early,caption="Optimal early exercise at time 2")
```

Table 4 implies that it is optimal to exercise the option at time t for the for all in-the-money paths. This leads to the following matrix that shows the cash flows conditional on not exercising prior to time 2. 

```{r echo=FALSE}
empty1 <-  rename(as.data.frame(matrix("--",8,1)),"t=1"="V1")
t3 <- rename(as.data.frame(matrix("0",8,1)),"t=3"="V1")
e2 <- opt_early$Exercise
cf2 <- cbind(empty1,rename(cbind(e2,t3),"t=2"="e2"))
cf2[3,2]=0
cf2[7,2]=0
kable(cf2, caption="Cash-flow matrix at time 2")
```

By recursively repeating the process for each time step I get the single optimal stopping time for each simulation path. It is then possible to average over all paths and compute the value of the American option. 



## Simulation of sample paths
An important first step in the Least Square Monte Carlo method is to generate a large number of stock price paths. It is reasonable to assume that the price paths follow a geometric Brownian motion REF Longstaff...  

The Brownian motion is a stochastic process $W(t)$ with the following properties:
\begin{enumerate}
\item $W(0)=0$
\item The mapping $t \rightarrow W(t)$ is with probability 1, a continous function on $[0,T]$
\item The increments $W(t_1)-W(t_0),W(t_2)-W(t_1),...,W(t_k)-W(t_{k-1})$ are independent for any k and any $0\leq t_0 < t_1 < ...< t_k \leq T$
\item $W(t)-W(s)\sim N(0,t-s)$ for any $0\leq s<t\leq T$
\end{enumerate}

The geometric Brownian motion is then simply a exponentiated Brownian motion with initial value $\log S(0)$. We specify the process of geometric Brownian motion as a stochastic derivative equation (SDE) in the form
$$\frac{dS(t)}{S(t)}=\mu \: d_t +\sigma \: dW(t),$$
where $\mu$ is the drift and $\sigma$ is the volatility of $\log S(t)$. It is not possible to write a solution $S(t)$ in terms of $W(t)$ but by applying Itô's lemma we can find the mean and higher order moments of the process. This gives us
$$d\: \log S(t)=(\mu-\frac{1}{2}\sigma^2)dt+\sigma \: dW(t)$$
where the $-\frac{1}{2}\sigma^2$ part corresponds to the difference between the median and the mean of the log-normal distribution. If $S(t)$ is a geometric Brownian motion distributed with mean $\mu$ and variance $\sigma^2$ and if S has initial value $S(0)$, then
$$S(t)=S(0)\exp \left([r-\frac{1}{2}\sigma^2]t+\sigma W(t) \right)$$
Given that the increments of W are independent and normally distributed we can write the following recursive procedure for generating price stock paths S, at $0=t_0<t_1<,...,t_n$:
$$S(t_{t+1})=S(t_i)\exp([r-\frac{1}{2}\sigma^2](t_{t+1}-t_i)0\sigma\sqrt{t_{i+1}-t_i}Z_{i+1})$$
Figure 1 below shows four simulations for different values of $\mu$ and $\sigma$. Having in mind that the owner of the put option only gains whenever the stock price is above the strike price, it seems that a large volatility is especially important for generating large payoffs. 

```{r, echo=FALSE, fig.cap="Simulation of stock price paths"}
gbm_vec <- function(nsim = 100, t = 25, mu = 0, sigma = 0.1, S0 = 100, dt = 1./365, type="normal") {
  
  # matrix of random draws - one for each day for each simulation
  if(type=="normal"){
    epsilon <- matrix(rnorm(t*nsim), ncol = nsim, nrow = t)
  }
  if(type=="sobol"){
    epsilon <- rnorm.sobol(t, nsim)
  }
  if(type=="halton"){
    epsilon <- rnorm.halton(t,nsim,init=FALSE)
  }
  
  
  # get GBM and convert to price paths
  gbm <- exp((mu - sigma * sigma / 2) * dt + sigma * epsilon * sqrt(dt))
  gbm <- apply(rbind(rep(S0, nsim), gbm), 2, cumprod)
  
  return(gbm)
}


nsim <- 25
t <- 100
mu <- 0.03
sigma <- 0.1
S0 <- 100

plot_gbm<- function(type="normal", mu=0.06, sigma=0.2){
  gbm <- gbm_vec(nsim, t, mu, sigma, S0, type=type)
  gbm_df <- as.data.frame(gbm) %>%
    mutate(ix = 1:nrow(gbm)) %>%
    pivot_longer(-ix, names_to = 'sim', values_to = 'price')
  gbm_df %>%
    ggplot(aes(x=ix, y=price, color=sim)) +
    geom_line() +
    theme(legend.position = 'none') +
    ggtitle(bquote(mu == .(mu) ~" and "~ sigma == .(sigma)))+
    ylim(55,150)+
    xlab("Time")+
    ylab("Price")
  
}
set.seed(3)
#plot_gbm("halton")
plot_1 <- plot_gbm("normal", mu=0.02, sigma=0.4)
plot_2 <- plot_gbm("normal", mu=0.08, sigma=0.4)
plot_3 <- plot_gbm("normal", mu=0.02, sigma=0.2)
plot_4 <- plot_gbm("normal", mu=0.08, sigma=0.2)
grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol=2)
```


# Variance reduction methods
In order to get an acceptable estimate of the option value a large number of simulations are needed. One possible solution is the use of variance reduction methods. 

## Antithetic variates
The method of antithetic variates introduces a perfectly negative correlation between pairs of simulation. 

$dS_{1,t} = rdS_{1,t}dt+\sigma dS_{1,t}dz_t$
$dS_{2,t} = rdS_{2,t}dt-\sigma dS_{2,t}dz_t$

This produces the pairs: $(Y_1,\tilde{Y}_1),(Y_2,\tilde{Y}_2),...,(Y_n\tilde{Y}_n)$ which are i.i.d. Furthermore, the pairs feature the same distribution, though they are not independent. 

Then, the antithetic variates estimator is simply the average of all 2n observations, 
\begin{equation}
  \label{eq:avv}
  \hat{Y}_{AV}=\frac{1}{2n}\left(\sum^n_{i=1}Y_i+\sum^n_{i=1}\tilde{Y}_i\right)=\frac{1}{n}\sum^n_{i=1}\left(\frac{Y_1+\tilde{Y}_1}{2}\right)
\end{equation}

It is apparent from the rightmost side of equation \ref{eq:avv} that $\hat{Y}_{AV}$ is the sample mean of the n independent observations 
\begin{equation}
\label{eq:obs}
  (\frac{Y_1+\tilde{Y}_1}{2}),(\frac{Y_2+\tilde{Y}_2}{2}),...,(\frac{Y_n+\tilde{Y}_n}{2})
\end{equation}

Thus, the central limit theorem can be applied and gives
$$\frac{\hat{Y}_{AV}-E[Y]}{\sigma_{AV}/\sqrt{n}}\Rightarrow N(0,1)$$
with 
$$\sigma^2_{AV}=Var \left[ \frac{Y_i+\tilde{Y_i}}{2} \right]$$

One might wonder when the antithetic variate estimator is to be preferred to an ordinary Monte Carlo estimator. We assume that the computational cost of computing a pair $(Y_i,\tilde{Y}_i)$ is approximately twice that of computing $Y_i$. If reducing variance is the main interest then using antithetics is preferred if
$$Var\left[\hat{Y}_{AV}\right]<Var\left[\frac{1}{2n}\sum^{2n}_{i=1}Y_i\right]$$
i.e., if
$$Var\left[Y_i+\tilde{Y_i} \right]<2Var[Y_i]$$
where we can use that $Y_i$ and $\tilde{Y_i}$ have the same variance if they have the same distribution and write the variance on the left as
$$Var[Y_i+\tilde{Y_i}]=Var[Y_i]+Var[\tilde{Y_i}]+2Cov[Y_i,\tilde{Y_i}]=2Var[Y_i]+2Cov[Y_i,\tilde{Y_i}]$$
Thus, the condition for antithetic variate to reduce variance becomes
$$Cov[Y_i,\tilde{Y_i}]<0$$
This requires that negative dependence in the inputs produce negative correlation between the outputs of paired replications. One example of this is when the option's payoff is a monotonous function of the underlying asset price e.g. an american put option that decreases in value whenever the value of the underlying increases. 



## Control variates
Control variates is yet another popular technique for reducing the variance of Monte Carlo simulation. The idea is to decompose the unknown expectation into the part known in closed form, and the part that needs to be estimated by simulation. Similar to the method of antithetic variates we generate another variable $X_i$ on each replication which produces a sequence of i.i.d. pairs $(X_i,Y_i)$ with a main objective of estimating the expected discounted payoff $E[Y_i]$. The usual estimator for $E[Y_i]$ is the sample mean $\bar{Y}=(Y_1+...+Yn)/n$, however there is an alternative. If we assume that the expectation E[X] of the $X_i$ is known then we can calculate, for any fixed coefficient b
$$Y_i(b)=Y_i-b(X_i-E[X])$$
$\lim_{n\rightarrow \infty}$
from the ith replication and the compute the sample mean
\begin{equation}
  \label{eq:cv_est}
  \bar{Y}(b)=\bar{Y}-b(\bar{X}-E[X])=\frac{1}{n}\sum^n_{i=1}(Y_i-b(X_i-E[X]))
\end{equation}
This is a control variate estimator because the observed error $\bar{X}-E[X]$ serves as a control in estimating $E[Y]$. Furthermore, we see that \ref{eqn:cv_est} is an unbiased estimator of $E[Y]$ because
$$E[\bar{Y}(b)]=E[\bar{Y}-b(\bar{X}-E[X])]=E[\bar{Y}]=E[Y]$$
and converges with probability 1 as $n\rightarrow \infty$,
\begin{equation}
\begin{split}
\lim_{n\rightarrow \infty}\frac{1}{n}\sum^n_{i=1}Y_i(b) & 
= \lim_{n\rightarrow \infty}\frac{1}{n}\sum^n_{i=1}Y_i(b) \\
& 
= E[Y-b(X-E[X])] \\
&
= E[Y]
\end{split}
\end{equation}

where the variance of $Y_i(b)$ is 

\begin{equation}
\begin{split}
Var[Y_i(b)]=Var[Y_i-b(X_i-E[X])]
& 
= \sigma^2_Y-2b\sigma_X\sigma_Y\rho{XY}+b^2\sigma_X^2 \equiv \sigma^2(b)
\end{split}
\end{equation}
The variance of the control variate estimator $\bar{Y}(b)$ is $\sigma^2(b)/n$ and the regular sample mean $\bar{Y}$ (which corresponds to b=0) has variance $\sigma^2_Y/n$. Hence, the control variate estimator has lower variance than the standard estimator if 
$b^2\sigma_X<2b\sigma_Y\rho_{XY}$


# Quasi-monte carlo methods
## Sobol draws
## Halton draws

# Results


