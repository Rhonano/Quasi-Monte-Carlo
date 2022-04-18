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
abstract: "Enter the text of your abstract here."
bibliography: citations.bib
biblio-style: unsrt
output: 
  pdf_document: 
    number_sections: true
    extra_dependencies: ["float"]
toc: true
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
```


```{r, echo=FALSE}
library(fOptions)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(LSMonteCarlo)
library(Jdmbs)
library(knitr)
library(bookdown)
library(gridExtra)
library(SobolSequence)
library(RQuantLib)
library(reshape2)
library(kableExtra)

 
QuasiAmerPutLSM <- function(Spot = 1, sigma = 0.2, n = 25000, m = 365, Strike = 1.1, 
                            r = 0.06, dr = 0, mT = 1, innovation="normal") 
{
  set.seed(25)
  if(innovation=="normal"){
    GBM <- rnorm.pseudo(n = n, dimension = m)
    for (i in 1:n) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    GBM[i, ])))}
  }
  if(innovation=="sobol"){
    GBM <- rnorm.sobol(n, m, scrambling = 1)
    for (i in 1:n) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    GBM[i, ])))}
  }
  if(innovation=="halton"){
    GBM <- rnorm.halton(n, m)
    for (i in 1:n) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    GBM[i, ])))}
  }
  
  X <- ifelse(GBM < Strike, GBM, NA)
  CFL <- matrix(pmax(0, Strike - GBM), nrow = n, ncol = m)
  Xsh <- X[, -m]
  X2sh <- Xsh * Xsh
  X3sh <- Xsh * Xsh * Xsh 
  Y1 <- CFL * exp(-1 * r * (mT/m))
  Y2 <- cbind((matrix(NA, nrow = n, ncol = m - 1)), Y1[, m])
  CV <- matrix(NA, nrow = n, ncol = m - 1)
  try(for (i in (m - 1):1) {
    reg1 <- lm(Y2[, i + 1] ~ Xsh[, i] + X2sh[, i]+ X3sh[, i])
    CV[, i] <- (matrix(reg1$coefficients)[1, 1]) + 
      ((matrix(reg1$coefficients)[2,1]) * Xsh[, i]) + 
      ((matrix(reg1$coefficients)[3,1]) * X2sh[, i]) +
    ((matrix(reg1$coefficients)[4,1]) * X3sh[, i])
    
    CV[, i] <- (ifelse(is.na(CV[, i]), 0, CV[, i]))
    Y2[, i] <- ifelse(CFL[, i] > CV[, i], Y1[, i], Y2[, i + 
                                                        1] * exp(-1 * r * (mT/m)))
  }, silent = TRUE)
  CV <- ifelse(is.na(CV), 0, CV)
  CVp <- cbind(CV, (matrix(0, nrow = n, ncol = 1)))
  POF <- ifelse(CVp > CFL, 0, CFL)
  FPOF <- firstValueRow(POF)
  dFPOF <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:m) {
    dFPOF[, i] <- FPOF[, i] * exp(-1 * mT/m * r * i)
  }
  PRICE <- mean(rowSums(dFPOF))
  SE <- sqrt(var(rowSums(dFPOF))/n)
  res <- list(price = (PRICE), Spot, Strike, sigma, n, m, r, 
              dr, mT)
  class(res) <- "AmerPut"
  return(c(PRICE,SE))
}



QuasiAmerPutLSM_AV <- function (Spot = 1, sigma = 0.2, n = 25000, m = 365, Strike = 1.1, 
                                r = 0.06, dr = 0, mT = 1, innovation="normal") 
{
  set.seed(25)
  if(innovation=="normal"){
    Z <- rnorm.pseudo(n = n, dimension = m)
    aZ <- Z * (-1)
    ZZ <- rbind(Z, aZ)
    GBM <- matrix(NA, nrow = n * 2, ncol = m)
    for (i in 1:(n * 2)) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    ZZ[i, ])))}
  }
  if(innovation=="sobol"){
    Z <- rnorm.sobol(n, m, scrambling = 1)
    aZ <- Z * (-1)
    ZZ <- rbind(Z, aZ)
    GBM <- matrix(NA, nrow = n * 2, ncol = m)
    for (i in 1:(n * 2)) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    ZZ[i, ])))}
  }
  if(innovation=="halton"){
    Z <- rnorm.halton(n, m)
    aZ <- Z * (-1)
    ZZ <- rbind(Z, aZ)
    GBM <- matrix(NA, nrow = n * 2, ncol = m)
    for (i in 1:(n * 2)) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    ZZ[i, ])))}
  }
  X <- ifelse(GBM < Strike, GBM, NA)
  CFL <- matrix(pmax(0, Strike - GBM), nrow = n * 2, ncol = m)
  Xsh <- X[, -m]
  X2sh <- Xsh * Xsh
  X3sh <- Xsh * Xsh * Xsh 
  Y1 <- CFL * exp(-1 * r * (mT/m))
  Y2 <- cbind((matrix(NA, nrow = n * 2, ncol = m - 1)), Y1[, m])
  CV <- matrix(NA, nrow = n * 2, ncol = m - 1)
  try(for (i in (m - 1):1) {
    reg1 <- lm(Y2[, i + 1] ~ Xsh[, i] + X2sh[, i]+X3sh[, i])
    CV[, i] <- (matrix(reg1$coefficients)[1, 1]) + 
      ((matrix(reg1$coefficients)[2,1]) * Xsh[, i]) + 
      ((matrix(reg1$coefficients)[3,1]) * X2sh[, i]) +
      ((matrix(reg1$coefficients)[4,1]) * X3sh[, i])
    CV[, i] <- (ifelse(is.na(CV[, i]), 0, CV[, i]))
    Y2[, i] <- ifelse(CFL[, i] > CV[, i], Y1[, i], Y2[, i + 
                                                        1] * exp(-1 * r * (mT/m)))
  }, silent = TRUE)
  CV <- ifelse(is.na(CV), 0, CV)
  CVp <- cbind(CV, (matrix(0, nrow = n * 2, ncol = 1)))
  POF <- ifelse(CVp > CFL, 0, CFL)
  FPOF <- firstValueRow(POF)
  dFPOF <- matrix(NA, nrow = n * 2, ncol = m)
  for (i in 1:m) {
    dFPOF[, i] <- FPOF[, i] * exp(-1 * mT/m * r * i)
  }
  PRICE <- mean(rowSums(dFPOF))
  table1 <- dFPOF[1:n,]
  table2 <- dFPOF[(n+1):(2*n),]
  meantable <- 0.5*(table1+table2)
  SE <- sqrt(var(rowSums(meantable))/n)
  res <- list(price = (PRICE), Spot, Strike, sigma, n, n, m, 
              r, dr, mT)
  class(res) <- "AmerPutAV"
  return(c(PRICE,SE))
}


QuasiAmerPutLSM_CV <- function (Spot = 1, sigma = 0.2, n = 25000, m = 365, Strike = 1.1, 
                                r = 0.06, dr = 0, mT = 1, innovation="normal") 
  
{
  set.seed(25)
  if(innovation=="normal"){
    GBM <- rnorm.pseudo(n = n, dimension = m)
    for (i in 1:n) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    GBM[i, ])))}
  }
  if(innovation=="sobol"){
    GBM <- rnorm.sobol(n, m, scrambling = 1)
    for (i in 1:n) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    GBM[i, ])))}
  }
  if(innovation=="halton"){
    GBM <- rnorm.halton(n, m)
    for (i in 1:n) {
      GBM[i, ] <- Spot * exp(cumsum(((r - dr) * (mT/m) - 0.5 * 
                                       sigma * sigma * (mT/m)) + (sigma * (sqrt(mT/m)) * 
                                                                    GBM[i, ])))}
  }
  
  
  X <- ifelse(GBM < Strike, GBM, NA)
  CFL <- matrix(pmax(0, Strike - GBM), nrow = n, ncol = m)
  Xsh <- X[, -m]
  X2sh <- Xsh * Xsh
  X3sh <- Xsh * Xsh * Xsh 
  Y1 <- CFL * exp(-1 * r * (mT/m))
  Y2 <- cbind((matrix(NA, nrow = n, ncol = m - 1)), Y1[, m])
  CV <- matrix(NA, nrow = n, ncol = m - 1)
  try(for (i in (m - 1):1) {
    reg1 <- lm(Y2[, i + 1] ~ Xsh[, i] + X2sh[, i]+X3sh[, i])
    CV[, i] <- (matrix(reg1$coefficients)[1, 1]) + 
      ((matrix(reg1$coefficients)[2,1]) * Xsh[, i]) + 
      ((matrix(reg1$coefficients)[3,1]) * X2sh[, i]) +
      ((matrix(reg1$coefficients)[4,1]) * X3sh[, i])
    
    CV[, i] <- (ifelse(is.na(CV[, i]), 0, CV[, i]))
    Y2[, i] <- ifelse(CFL[, i] > CV[, i], Y1[, i], Y2[, i + 
                                                        1] * exp(-1 * r * (mT/m)))
  }, silent = TRUE)
  CV <- ifelse(is.na(CV), 0, CV)
  CVp <- cbind(CV, (matrix(0, nrow = n, ncol = 1)))
  POF <- ifelse(CVp > CFL, 0, CFL)
  FPOF <- firstValueRow(POF)
  dFPOF <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:m) {
    dFPOF[, i] <- FPOF[, i] * exp(-1 * mT/m * r * i)
  }
  PRICE1 <- mean(rowSums(dFPOF))
  dCFL <- CFL[, m] * exp(-1 * r * mT)
  EUsimul <- mean(dCFL)
  EUBS <- EuPutBS(Spot, sigma, Strike, r, dr, mT)
  PRICE <- PRICE1 - (EUsimul - EUBS)
  SE <- sqrt(var(rowSums(dFPOF) - (dCFL - EUBS))/n)
  res <- list(price = (PRICE), Spot, Strike, sigma, n, m, r, 
              dr, mT)
  class(res) <- "AmerPutCV"
  return(c(PRICE,SE))
}







```

\newpage

# Introduction

One of the most important option pricing formula is the @blackscholes which offers a closed-form solution to European options where only the asset price is a stochastic factor. Since then more flexibility has been added via finite difference and binomial models such as the one proposed by @cox1979. However, when contingencies become more complex, such as the possibility of early exercise, there exist no analytical solution, and this poses a problem as many of the most traded options are American styled^[https://www.tradestation.com/insights/top-stocks-for-options-trading/]. Fortunately, Monte Carlo simulation can be used to represent the randomness of stock price paths and is readily applied when multiple factors are stochastic. Specifically, the least squares Monte Carlo (LSM) approach by @longstaff estimates the conditional expected cash flow from a least squares regression using simulated cross-sectional price paths. The standard LSM algorithm has some flaws such as requiring many observations and thus a large computation time in order to converge to the true price. This paper attempts to price Bermudan-style american put options and improve the LSM by regressing a cubic basis function and applying variance reduction methods in the form of antithetic and control variates as proposed by @hammersley1956 and @hull_white_1988. Furthermore, we compare standard pseudo-random numbers to quasi-random (low-discrepancy) Sobol numbers developed by @sobol1967 in the hopes that it will reduce the asymptotic error and oscillations of the price estimates.  

# American option valuation problem
In order to understand how to valuate American put options it is necessary to define how they differ from the more simple European options. A put option is a security that gives the owner the right to sell a fixed number of shares of an underlying stock at a fixed price. There are many types of options where the most common two are European and American. The European type can only be exercised on the maturity date whereas the American type can be exercised at any time up until the expiration. This extra flexibility means that the price of an American option will always be equal to or larger than that of the European. In this paper we simplify the computation process by pricing Bermudan-style American options which can only be exercised once a day. Furthermore, we distinguish between options that are out-of-the-money (OTM), at-the-money (ATM) and in-the-money (ITM). OTM means that the underlying asset price is below the strike price e.g., the intrinsic value $(S_T-K)_+$ is zero (the time value that we attempt to price, however, is usually still larger than 0). ATM and ITM options have strike prices that are at and above the underlying asset price, respectively. ITM options should always be more expensive than ATM and OTM due to the larger intrinsic value, however, there might be irregularities due to poor liquidity and large bid-ask spreads. The pricing problem of American options can be simplified to specifying a process $U(t), 0\leq t \leq T$ which represents the discounted payoff from exercise at time t, and a class of stopping times $\mathcal T$. The problem is then to find the optimal expected discounted payoff$$\sup_{\tau\in \mathcal T}E[U(\tau)]$$

# Data
No data is needed when option pricing with LSM Monte Carlo because it is based on simulations. However, it can be useful as a benchmark when comparing performance of different models. There are several ways to choose a benchmark for your results like comparing it to the finite differences or closed form European price as done by @longstaff. We choose to compare the results against observed option prices. In order to get consistent results it is important to compare with liquid options that do not have too large bid-ask spreads. For that reason we choose a 10 week maturity, slightly ITM SPDR S&P 500 ETF (SPY) put for the numerical calculations. Another important input for the Monte Carlo simulation is to choose the right volatility. One option is to base the volatility on a historical mean, however, the implied volatility usually ranges between $10-30\%$ and can go as high as $58\%$ as it did during the Covid-19 crisis^[https://www.alphaquery.com/stock/SPY/volatility-option-statistics/30-day/iv-mean] and choosing the wrong volatility will affect the price (see figure 1). We use the alternative of calculating the option price with finite differences method and calculating the implied volatility numerically. This method gives a great estimate on what volatility the true price is based on. 
<br/>
<br/>
Furthermore, for discounting the cash flows we assume that the risk-free rate is $2.72\%$^[https://ycharts.com/indicators/10_year_treasury_rate] and that SPY pays $1.228\%$ in dividends per year as of 12/04/2022^[https://www.macrotrends.net/stocks/charts/SPY/spdr-s-p-500-etf/dividend-yield-history]. Table 1 shows the observed and estimated prices of a slightly ITM SPY PUT (strike price of 450 with an asset price of 447.57). The relationship between price, strike and volatility is visualized with a surface plot figure 1 where it is clear to see that both volatility and strike price have significant positive effects on the price of the option. 
```{r, echo=FALSE}
price_10 <- QuasiAmerPutLSM(Spot = 447.57, sigma = 0.10, n = 50000, m = 70, Strike = 450,
                     r = 0.0227, dr = 0.01228, mT = 70/365, innovation="normal")[1]
price_30 <- QuasiAmerPutLSM(Spot = 447.57, sigma = 0.30, n = 50000, m = 70, Strike = 450,
                     r = 0.0227, dr = 0.01228, mT = 70/365, innovation="normal")[1]
true_price <- 16.11
df<- data.frame(price_10,true_price,price_30)
colnames(df) <- c("10% volatility", "Observed price", "30% volatility")
kable(df, digits=2, caption="Effect of volatility on option price")
```

```{r echo=FALSE, fig.height=8, fig.width=8, cache=TRUE}
surface <- AmerPutLSMPriceSurf(Spot=447.57, vols=seq(0.1,0.4,0.03), n=10000, m=70, strikes=seq(420,475,4), r=0.0227, dr=0.01228, mT=70/365)
plot(surface, color = divPalette(150, "RdBu"), caption="Figure 1: 3D plot of option price surface")
```


# Least squares Monte Carlo

## Least Square Monte Carlo algorithm
The LSM approach uses least squares regression to recursively approximate the discounted payoff at $t_{k-1},t_{k-2},...,t_1$. The reason for working backwards is that the cash flows can differ since it may be optimal to stop before the maturity date. We assume that the unknown functional form of $F(\omega,t_{K-1})$ for approximating the conditional expectation can be represented as a linear combination of basis functions. Assuming that the value of the asset underlying the option follows a Markov process^[A Markov chain is a stochastic model describing a sequence of possible events in which the probability of each event depends only on the state attained in the previous event @gagniuc2017]. We let the choice of basis function be the cubic polynomial, such that $F(\omega,t_{K-1})=a_0X+a_1X^2+a_2X^3$, where $a_0$, $a_1$ and $a_2$ are constant coefficients. To implement the LSM algorithm, we then regress the discounted payoffs onto the basis functions where the option is in the money at time $t_{k-1}$. The advantage of only regressing on paths that are in the money is that we can limit the region over which the conditional expectation must be estimated. Since the values of the basis function are i.i.d., we can use theorem 3.5 of @white2014 to show that fitted value of the regression $\hat{F}(\omega,t_{K-1})$ converge in probability to $F(\omega,t_{K-1})$ as the number of in-the-money price paths, N, goes to infinity. Furthermore, Theorem 1.2.1 of @takeshi1985 implies that $\hat{F}(\omega,t_{K-1})$ is the best linear unbiased estimator of $F(\omega,t_{K-1})$. 
<br/>
<br/>
Once the conditional expectation function at time $t_{k-1}$ is estimated we can determine whether it is optimal to exercise early for path $\omega$ by comparing the immediate exercise value to $\hat{F}(\omega,t_{K-1})$. Once the optimal exercise time is decided the option cash flow can be estimated. The option value is then found by starting at time zero, moving forward along each path $\omega$ until each optimal stopping time has been reached. We then discount the cash flow from exercise back to time zero and averaging over all paths. For a simple numerical example of the LSM algorithm see the appendix. 


## Simulation of sample paths
An important first step in the Least Square Monte Carlo method is to generate a large number of stock price paths. It is reasonable to assume that the price paths follow a geometric Brownian motion as in @longstaff.   

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
Figure 2 below shows four scenarios for different values of $\mu$ and $\sigma$. Having in mind that the owner of the put option only gains whenever the stock price is below the strike price (and has limited loss when above), it seems that a large volatility is especially important for generating large payoffs while a low risk-free rate helps drive the overall trend down (or less up). 

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
grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol=2, top="Figure 2: Price path simulations")
```


# Variance reduction methods
In order to get an acceptable estimate of the option value a large number of simulations are needed. One possible solution is the use of variance reduction methods. In this paper we use the antithetic and control variates in context of quasi-Monte Carlo methods. 

## Antithetic variates
The method of antithetic variates introduces a perfectly negative correlation between pairs of simulation. 

$$dS_{1,t} = rdS_{1,t}dt+\sigma dS_{1,t}dz_t$$
$$dS_{2,t} = rdS_{2,t}dt-\sigma dS_{2,t}dz_t$$

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
This is a control variate estimator because the observed error $\bar{X}-E[X]$ serves as a control in estimating $E[Y]$. Furthermore, we see that equation 3 is an unbiased estimator of $E[Y]$ because
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


# Quasi-Monte Carlo methods
A requirement for using Monte-Carlo simulation is to sample from a suitable random distribution. However, computers are not capable of generating truly random numbers. Instead they generate pseudo-random numbers which, when enough are considered together, appear random. Relying on pseudo-random generated numbers and using statistical properties of assumed random sequences has two drawbacks. First, the validity of pseudo-random numbers is somewhat questionable. Second, the rate of convergence decreases slowly and more draws are needed. 

## Sobol draws
One popular method of generating quasi-random numbers is the Sobol sequence, developed by @sobol1967 and was improved by @antonov1979. The Sobol sequence is a collection of D sequences of numbers in (0,1). The sequences consist of integers $(y_i^d)$ between 0 and $2^{32}-1$, so the $i$th number of the axis $d$ is $$x_i^d=\frac{y^d_i}{2^{32}}\in(0,1)$$

The integers $y^d_i$ on a given axis d are produced with $y_0^d=0$ and
$$y^d_{i+1}=y_i^d\oplus DN^d_{J_i}$$
Where $\oplus$ denotes the bitwise exclusive or operator, $J_i$ is the rightmost 0 bit in the binary expansion of $i$ and $\{DN_j^d,0\leq j\leq31\}$ are freely chosen numbers that initialize the sequence. The numbers are generated in the same way as in @savine2018

\begin{equation}
\begin{split}
x_0^0 & = 0 \\
x_1^0 & = (0.1)_2 = 0.5 \\
x_2^0 & = (0.11)_2 = 0.75 \\
x_3^0 & = (0.1)_2 \oplus  (0.11)_2 = (0.01)_2= 0.25 \\
x_4^0 & = (0.111)_2 = 0.875 \\
x_5^0 & = (0.1)_2 \oplus  (0.111)_2 = (0.011)_2= 0.375
\end{split}
\end{equation}

Thus, the Sobol sequence increases uniformity by generating numbers that fill out the spaces with the largest gap. This is also apparent when plotting the pseudo random numbers and Sobol sequence in two dimensions

```{r echo=FALSE}

srange <- sobolSequence.dimMinMax()
mrange <- sobolSequence.dimF2MinMax(srange[1])
points <- sobolSequence.points(dimR=srange[1], dimF2=mrange[1], count=512)

sobol_plot <- as.data.frame(points) %>%
  ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  ggtitle("2D Sobol", subtitle="sample size=512")

pseudo_plot <- as.data.frame(runif.pseudo(n = 512, dimension = 2)) %>% 
  ggplot(aes(x=V1, y=V2))+
  geom_point()+
  ggtitle("2D Pseudo", subtitle="sample size=512")

grid.arrange(sobol_plot, pseudo_plot, ncol=2, top="Figure 3: Comparison of random draws")

```

The final step before we can use the sequence in the Monte Carlo simulation is to transform the uniform distribution to a Gaussian. This is done using the Box-Muller transform which takes two independent uniform samples and sends them to two independent random variables with Gaussian distribution. For more information see @blitzstein2015

The main benefit of using the Sobol sequence is that it offers a lower discrepancy, meaning that the space of possibilities are spread more evenly. This generally results in a faster convergence speed of $O\left( \frac{1}{N}\right)$, rather than $O \left( \frac{1}{\sqrt{N}} \right)$ gained from pseudo-random sequences. The results are shown in the following chapter. 

# Numerical results
First, it can be interesting to take a look at how the different variance reducing model are pricing options compared to the observed data. Table 2 and 3 are shown for SPY put options with a 70 day maturity, an underlying asset price of 445.57, a dividend yield of 1.28%, a risk-free rate of 2.27% and 50,000 simulations. In table 2 we see no large differences in price estimates due to the large number simulations, however, it still seems that the Sobol control and antithetic variates are pricing closer to the observed than the rest. Comparing the different estimates to the observed price with such as small sample should only be taken with a grain of salt as the observed prices can differ a lot due to varying volumes and bid-ask spread sizes. However, the errors across all models range between 0-10 cents (0%-0.71% relative error) which is not bad considering that the bid-ask spread is sometimes several dollars (or 5-10% of the price). It can also be noted that in 8/9 of the estimates the Sobol version prices closer to the observed price than the pseudo version. 
<br/>
<br/>
Table 3 better shows how control and antithetic variates can be used to decrease standard errors of the estimator. Interestingly, we see that the control variates perform better for OTM options (low strike prices) whereas the antithetic variates help decrease the standard error for ITM options. This is also consistent with the results of @boire2021. The curious behavior of the control variates is discussed in @moreni2003. In short, the explanation is that as the difference between the prices of American and European options (the exercise premium) gets larger, the change of measure that would be optimal for the European option gets further away from the one that would be optimal for the American option. 


```{r echo=FALSE, cache=TRUE}
# 447.57 closing price per 08/04/2022, exp 2022-06-17. mat=70, libor/rf=2.27


otm_vol <- AmericanOptionImpliedVolatility(type="put", value=9.92, underlying=447.57,
                                strike=430, dividendYield=0.0128, riskFreeRate=0.0227,
                                maturity=70/365, volatility=0.19)
# OTM Strike 430, mean(bid,ask)=9.92, iv=0.2278214
atm_vol <-AmericanOptionImpliedVolatility(type="put", value=16.11, underlying=447.57,
                                strike=450, dividendYield=0.0128, riskFreeRate=0.0227,
                                maturity=70/365, volatility=0.19)
# ATM Strike 450 mean(bid,ask)=16.11, iv=0.1918263
itm_vol <-AmericanOptionImpliedVolatility(type="put", value=26.95, underlying=447.57,
                                strike=470, dividendYield=0.0128, riskFreeRate=0.0227,
                                maturity=70/365, volatility=0.19)
# ITM Strike 470 mean(bid,ask)=26.95, iv=0.1673667

strikes<- c(430,450,470)
ivs<- c(otm_vol,atm_vol,itm_vol)
mat=70

norm_list <- list()
norm_list_SE <- list()
for(i in 1:length(strikes)){
  output <- QuasiAmerPutLSM(Spot = 447.57, sigma = ivs[i], n = 50000, m = mat, Strike = strikes[i], 
                            r = 0.0227, dr = 0.0128, mT = mat/365, innovation="normal")
  price <- output[1]
  se <- output[2]
  norm_list <- append(norm_list,price) 
  norm_list_SE <- append(norm_list_SE,se) 
}

sobol_list <- list()
sobol_list_SE <- list()
for(i in 1:length(strikes)){
  output <- QuasiAmerPutLSM(Spot = 447.57, sigma = ivs[i], n = 50000, m = mat, Strike = strikes[i], 
                            r = 0.0227, dr = 0.0128, mT = mat/365, innovation="sobol")
  price <- output[1]
  se <- output[2]
  sobol_list <- append(sobol_list,price) 
  sobol_list_SE <- append(sobol_list_SE,se) 
}

norm_cv_list <- list()
norm_cv_list_SE <- list()
for(i in 1:length(strikes)){
  output <- QuasiAmerPutLSM_CV(Spot = 447.57, sigma = ivs[i], n = 50000, m = mat, Strike = strikes[i], 
                            r = 0.0227, dr = 0.0128, mT = mat/365, innovation="normal")
  price <- output[1]
  se <- output[2]
  norm_cv_list <- append(norm_cv_list,price) 
  norm_cv_list_SE <- append(norm_cv_list_SE,se) 
}


sobol_cv_list <- list()
sobol_cv_list_SE <- list()
for(i in 1:length(strikes)){
  output <- QuasiAmerPutLSM_CV(Spot = 447.57, sigma = ivs[i], n = 50000, m = mat, Strike = strikes[i], 
                            r = 0.0227, dr = 0.0128, mT = mat/365, innovation="sobol")
  price <- output[1]
  se <- output[2]
  sobol_cv_list <- append(sobol_cv_list,price) 
  sobol_cv_list_SE <- append(sobol_cv_list_SE,se) 
}

norm_av_list <- list()
norm_av_list_SE <- list()
for(i in 1:length(strikes)){
  output <- QuasiAmerPutLSM_AV(Spot = 447.57, sigma = ivs[i], n = 50000, m = mat, Strike = strikes[i], 
                               r = 0.0227, dr = 0.0128, mT = mat/365, innovation="normal")
  price <- output[1]
  se <- output[2]
  norm_av_list <- append(norm_av_list,price) 
  norm_av_list_SE <- append(norm_av_list_SE,se) 
}

sobol_av_list <- list()
sobol_av_list_SE <- list()
for(i in 1:length(strikes)){
  output <- QuasiAmerPutLSM_AV(Spot = 447.57, sigma = ivs[i], n = 50000, m = mat, Strike = strikes[i], 
                               r = 0.0227, dr = 0.0128, mT = mat/365, innovation="sobol")
  price <- output[1]
  se <- output[2]
  sobol_av_list <- append(sobol_av_list,price) 
  sobol_av_list_SE <- append(sobol_av_list_SE,se) 
}



results_df <- data.frame(cbind("Strike Price"=strikes,
                               "Gaussian"=unlist(norm_list),
                               "Sobol"=unlist(sobol_list),
                               "Gaussian CV"=unlist(norm_cv_list),
                               "Sobol CV"=unlist(sobol_cv_list),
                               "Gaussian AV"=unlist(norm_av_list),
                               "Sobol AV"=unlist(sobol_av_list),
                               c(9.92,16.11,26.95))) %>% 
  round(3)


colnames(results_df) <- c("Strike price","Pseudo","Sobol","Pseudo control variates","Sobol control variates", "Pseudo antithetic variates","Sobol antithetic variates", "Observed")


results_df_SE <- data.frame(cbind("Strike Price"=strikes,
                               "Gaussian"=unlist(norm_list_SE),
                               "Sobol"=unlist(sobol_list_SE),
                               "Gaussian CV"=unlist(norm_cv_list_SE),
                               "Sobol CV"=unlist(sobol_cv_list_SE),
                               "Gaussian AV"=unlist(norm_av_list_SE),
                               "Sobol AV"=unlist(sobol_av_list_SE)
                               )) %>% 
  round(3)
colnames(results_df_SE) <- c("Strike price","Pseudo","Sobol","Pseudo control variates","Sobol control variates", "Pseudo antithetic variates","Sobol antithetic variates")


```

```{r, echo=FALSE}
kable(results_df, caption="Price estimates", digits=2) %>%
  kable_styling(latex_options = c("striped", "scale_down"))
  
kable(results_df_SE, caption="Standard errors", digits=3) %>%
  kable_styling(latex_options = c("striped", "scale_down"))
```


Figure 4 shows how the standard errors of the different models converge across the number of simulations for the ATM put option. 4.a shows the same story as table 3 - that is, antithetic variates perform best in regards to starndard errors with the control variates coming in at a close second. However, there seems to be differences between pseudo and Sobol draws that become more apparent when zooming in on the plot as in 4.b. The benefit of using Quasi-random draws such as Sobol is that the standard error is generally lower but also much more consistent. Looking at figure 5.b we also see that the Sobol needs much fewer simulations to converge to the true price (as indicated by the black dotted line) whereas the pseudo draws are still inconsistent, even at 50,000 simulations. 

```{r, echo=FALSE, cache=FALSE, cache=TRUE}
seq_list_1 <- seq(10,5000,500)
seq_list_2 <- seq(4000,21000,1000)
seq_list_3 <- seq(25000,50000,1500)
seq_list <- append(append(seq_list_1,seq_list_2), seq_list_3)




# 447.57 closing price per 08/04/2022, exp 2022-06-17. mat=70, libor/rf=2.27


otm_vol <- AmericanOptionImpliedVolatility(type="put", value=9.92, underlying=447.57,
                                           strike=430, dividendYield=0.0128, riskFreeRate=0.0227,
                                           maturity=70/365, volatility=0.19)
# OTM Strike 430, mean(bid,ask)=9.92, iv=0.2278214
atm_vol <-AmericanOptionImpliedVolatility(type="put", value=16.11, underlying=447.57,
                                          strike=450, dividendYield=0.0128, riskFreeRate=0.0227,
                                          maturity=70/365, volatility=0.19)
# ATM Strike 450 mean(bid,ask)=16.11, iv=0.1918263
itm_vol <-AmericanOptionImpliedVolatility(type="put", value=26.95, underlying=447.57,
                                          strike=470, dividendYield=0.0128, riskFreeRate=0.0227,
                                          maturity=70/365, volatility=0.19)
# ITM Strike 470 mean(bid,ask)=26.95, iv=0.1673667

strikes<- c(430,450,470)
ivs<- c(otm_vol,atm_vol,itm_vol)
mat=70

norm_list_c <- list()
norm_list_c_SE <- list()
for(i in seq_list){
  output <- QuasiAmerPutLSM(Spot = 447.57, sigma = atm_vol, n = i, m = mat, Strike = 450, 
                            r = 0.0227, dr = 0.0128, mT = mat/365, innovation="normal")
  price <- output[1]
  se <- output[2]
  norm_list_c <- append(norm_list_c,price) 
  norm_list_c_SE <- append(norm_list_c_SE,se) 
}

sobol_list_c <- list()
sobol_list_c_SE <- list()
for(i in seq_list){
  output <- QuasiAmerPutLSM(Spot = 447.57, sigma = atm_vol, n = i, m = mat, Strike = 450, 
                            r = 0.0227, dr = 0.0128, mT = mat/365, innovation="sobol")
  price <- output[1]
  se <- output[2]
  sobol_list_c <- append(sobol_list_c,price) 
  sobol_list_c_SE <- append(sobol_list_c_SE,se) 
}

norm_cv_list_c <- list()
norm_cv_list_c_SE <- list()
for(i in seq_list){
  output <- QuasiAmerPutLSM_CV(Spot = 447.57, sigma = atm_vol, n = i, m = mat, Strike = 450, 
                               r = 0.0227, dr = 0.0128, mT = mat/365, innovation="normal")
  price <- output[1]
  se <- output[2]
  norm_cv_list_c <- append(norm_cv_list_c,price) 
  norm_cv_list_c_SE <- append(norm_cv_list_c_SE,se) 
}


sobol_cv_list_c <- list()
sobol_cv_list_c_SE <- list()
for(i in seq_list){
  output <- QuasiAmerPutLSM_CV(Spot = 447.57, sigma = atm_vol, n = i, m = mat, Strike = 450, 
                               r = 0.0227, dr = 0.0128, mT = mat/365, innovation="sobol")
  price <- output[1]
  se <- output[2]
  sobol_cv_list_c <- append(sobol_cv_list_c,price) 
  sobol_cv_list_c_SE <- append(sobol_cv_list_c_SE,se) 
}

norm_av_list_c <- list()
norm_av_list_c_SE <- list()
for(i in seq_list){
  output <- QuasiAmerPutLSM_AV(Spot = 447.57, sigma = atm_vol, n = i, m = mat, Strike = 450, 
                               r = 0.0227, dr = 0.0128, mT = mat/365, innovation="normal")
  price <- output[1]
  se <- output[2]
  norm_av_list_c <- append(norm_av_list_c,price) 
  norm_av_list_c_SE <- append(norm_av_list_c_SE,se) 
}

sobol_av_list_c <- list()
sobol_av_list_c_SE <- list()
for(i in seq_list){
  output <- QuasiAmerPutLSM_AV(Spot = 447.57, sigma = atm_vol, n = i, m = mat, Strike = 450, 
                               r = 0.0227, dr = 0.0128, mT = mat/365, innovation="sobol")
  price <- output[1]
  se <- output[2]
  sobol_av_list_c <- append(sobol_av_list_c,price) 
  sobol_av_list_c_SE <- append(sobol_av_list_c_SE,se) 
}




results_df_c <- data.frame(cbind("Simulations"=seq_list,
                               "Gaussian"=unlist(norm_list_c),
                               "Sobol"=unlist(sobol_list_c),
                               "Gaussian CV"=unlist(norm_cv_list_c),
                               "Sobol CV"=unlist(sobol_cv_list_c),
                               "Gaussian AV"=unlist(norm_av_list_c),
                               "Sobol AV"=unlist(sobol_av_list_c),
                               16.11
                               )) 


colnames(results_df_c) <- c("Simulations","Pseudo","Sobol","Pseudo CV","Sobol CV", "Pseudo AV","Sobol AV", "Observed")


results_df_SE_c <- data.frame(cbind("Simulations"=seq_list,
                                  "Gaussian"=unlist(norm_list_c_SE),
                                  "Sobol"=unlist(sobol_list_c_SE),
                                  "Gaussian CV"=unlist(norm_cv_list_c_SE),
                                  "Sobol CV"=unlist(sobol_cv_list_c_SE),
                                  "Gaussian AV"=unlist(norm_av_list_c_SE),
                                  "Sobol AV"=unlist(sobol_av_list_c_SE)
)) 


colnames(results_df_SE_c) <- c("Simulations","Pseudo","Sobol","Pseudo CV","Sobol CV", "Pseudo AV","Sobol AV")


df.melted_se<-melt(results_df_SE_c, id="Simulations")
high_se <- ggplot(df.melted_se, aes(x=Simulations, y=value, color=variable))+
         geom_line()+
  ylim(0,0.4)+
  ggtitle("4.a: High number of simulations")+
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())+
  ylab("Standard deviation")+
  guides(color=guide_legend(nrow=3, byrow=TRUE))




norm_cv_list_c_2 <- list()
norm_cv_list_c_SE_2 <- list()
for(i in seq(100,2000,20)){
  output <- QuasiAmerPutLSM_AV(Spot = 447.57, sigma = atm_vol, n = i, m = mat, Strike = 450, 
                               r = 0.0227, dr = 0.0128, mT = mat/365, innovation="normal")
  price <- output[1]
  se <- output[2]
  norm_cv_list_c_2 <- append(norm_cv_list_c_2,price) 
  norm_cv_list_c_SE_2 <- append(norm_cv_list_c_SE_2,se) 
}


sobol_cv_list_c_2 <- list()
sobol_cv_list_c_SE_2 <- list()
for(i in seq(100,2000,20)){
  output <- QuasiAmerPutLSM_AV(Spot = 447.57, sigma = atm_vol, n = i, m = mat, Strike = 450, 
                               r = 0.0227, dr = 0.0128, mT = mat/365, innovation="sobol")
  price <- output[1]
  se <- output[2]
  sobol_cv_list_c_2 <- append(sobol_cv_list_c_2,price) 
  sobol_cv_list_c_SE_2 <- append(sobol_cv_list_c_SE_2,se) 
}

results_df_SE_c_2 <- data.frame(cbind("Simulations"=seq(100,2000,20),
                                    "Gaussian AV"=unlist(norm_cv_list_c_SE_2),
                                    "Sobol AV"=unlist(sobol_cv_list_c_SE_2)
                                    
)) 
results_df_c_2 <- data.frame(cbind("Simulations"=seq(100,2000,20),
                                      "Gaussian AV"=unlist(norm_cv_list_c_2),
                                      "Sobol AV"=unlist(sobol_cv_list_c_2)
                                      
)) 

colnames(results_df_SE_c_2) <- c("Simulations","Pseudo CV","Sobol CV")
colnames(results_df_c_2) <- c("Simulations","Pseudo CV","Sobol CV")
df.melted_2_SE<-melt(results_df_SE_c_2, id="Simulations")
df.melted_2<-melt(results_df_c_2, id="Simulations")

low_se <- ggplot(df.melted_2_SE, aes(x=Simulations, y=value, color=variable))+
  geom_line()+
  ggtitle("4.b: Low number of simulations")+
  theme(legend.position="bottom")+
  ylab("Standard deviation")

grid.arrange(high_se, low_se, ncol=2, top="Figure 4: Convergence of standard deviation")

```

Finally, looking at figure 5.a we see that even though the control variate method has lower standard deviation than the standard method, it still requires just as many (around 25k) simulations in order to converge. From the colored shades that represent the standard error bands, we observe that both the control and antithetic variates are more precise than the standard method. Between figure 4 and 5 it is clear to see that the quasi-random Monte Carlo methods with antithetic variates perform best both in regards to standard deviation, consistency and speed of convergence. These are similar findings as to @birge1995 who find that Sobol draws consistently require fewer observations to converge and have less fluctuations than the standard Monte Carlo. The only disadvantage is that is has a larger random number generation time, however, that is still outweighed by a much lower total computation time. 

```{r, echo=FALSE, cache=FALSE, cache=TRUE}
df.melted<-melt(results_df_c, id="Simulations")
total.melted <- merge(df.melted,df.melted_se,by=c("Simulations","variable"))

total.melted_filter <- total.melted[total.melted$variable=="Pseudo" | total.melted$variable=="Sobol",]
quasiplot <- ggplot(total.melted_filter, aes(x=Simulations, y=value.x, color=variable, ymin=value.x-value.y, ymax=value.x+value.y, fill=variable))+
  geom_line()+
  ylim(15.5,17)+
  geom_ribbon( alpha=0.3, colour=NA)+
  theme(legend.position="bottom")+
  ylab("Price, USD")+
  geom_hline(yintercept=16.11, linetype="dashed", color = "black")+
  ggtitle("5.b: Quasi-random methods")


total.melted_filter_2 <- total.melted[total.melted$variable=="Sobol" | total.melted$variable=="Sobol AV" | total.melted$variable=="Sobol CV",]
varplot <- ggplot(total.melted_filter_2, aes(x=Simulations, y=value.x, color=variable, ymin=value.x-value.y, ymax=value.x+value.y, fill=variable))+
  geom_line()+
  ylim(15.5,17)+
  geom_ribbon( alpha=0.2,colour = NA)+
  theme(legend.position="bottom")+
  ylab("Price, USD")+
  geom_hline(yintercept=16.11, linetype="dashed", color = "black")+
  ggtitle("5.a: Variance reduction methods")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

grid.arrange(varplot, quasiplot, ncol=2, top="Figure 5: Convergence of price")
```



# Discussion and extensions
The LSM algorithm is flexible in terms of modifying either the GBM process, the least squares regression or adding variance reduction methods. The authors of the LSM paper @longstaff argue that you can use other polynomials than the simple cubic as basis functions for the regression. Longstaff and Schwartz use the weighted Laguerre polynomials $L_n(X)=\exp(-X/2)\frac{e^X}{n!}\frac{d^n}{dX^n}(X^ne^{-X})$ but find that Fourier and trigonometric series or even the simple quadratic polynomial can yield accurate results. @moreno2003 test for the robustness and convergence of different polynomials for the LSM algorithm and find that the technique is very robust relative to the type and number of basis functions. Choosing different polynomials yield similar results which is not surprising as we can express any polynomial as a linear combination of others. Furthermore they find that in- and out-of-sample option prices are alike and standard errors are low in both cases. Finally, using a reasonable amount of basis functions (between 3 and 20) also produce similar results and using more can even lead to numerical problems in the regression. 

Another problems lies not in the LSM algorithm but in the assumptions of the stock price path. Besides the continuous log-normally distributed component the geometric Brownian motion can be extended with a compound Poisson jump process to capture negative skewness and excess kurtosis of the stock price movements. This was first introduced with the Merton log-normal jump diffusion process @merton1976 and later @kou2002 with double exponential distributed jumps. @kou2004 use the double exponential diffusion model and compute maximum relative errors that are in most cases below $1\%$. This is in the range as the model presented in this paper so the jump extension boils down to choosing between analytical tractability and reality. They also argue that the jump diffusion process might be more suitable for the currency market than the equity market. 

# Conclusion
The standard Monte Carlo requires a large amount of simulations in order to converge to the true price and even then it is not stable. There are ways to improve the method ranging from variance reduction techniques to low discrepancy sequences and thankfully these are non-exclusive and can be used together to improve the model further. This paper finds that a combination of antithetic variates and quasi-random sobol draws decrease the standard error for all measured strike prices, increases the speed of convergence and reduces the number of fluctuations compared to the standard Monte Carlo. The control variate method can also reduce the standard error efficiently for OTM option, however, the effect declines as the strike price increases.  


# References

<div id="refs"></div>

\newpage
# Appendix

## A numerical example
The key idea of the LSM approach is to estimate the conditional expected discounted payoff from the cross-sectional information of the Monte Carlo simulation. 

For a better understanding we produce a simple numerical example below. Consider an American put option on a share of a non-dividend-paying stock. The option is exercisable at times $t=1,2,3$ where 3 is the maturity date. The paths are generated using a geometric Brownian motion with a risk free rate of 6%. 

```{r, echo=FALSE, cache=TRUE}
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
The objective is then to recursively compute intermediate matrices to solve for the stopping rule that maximizes the value of the option at each time t, beginning at $t=3$. These cash flows are also equal to that of the European option. 


```{r echo=FALSE, cache=TRUE}
empty <- rename(as.data.frame(matrix("--",8,2)),"t=1"="V1","t=2"="V2")
cflow <- pmax(0,1.05-num_ex["t=3"])
cflow3 <- cbind(empty, cflow)
cflow3[is.na(cflow3)]=0
cflow3 <- rename(cflow3,"t=3"="cflow")
kable(cflow3, caption="Cash flow at t=3")
```
The 0 values indicate that the path is out-of-money and not worth anything. By only using the paths that are in-the-money we get a better estimate of the conditional expectation function in addition to speeding up the algorithm. Let X denote the stock prices as time 2 and Y denote the discounted cash flows at time 3 if not exercised at time 2. 

```{r echo=FALSE, cache=TRUE}
Y <- round(cflow3["t=3"]/(1+(0.1)/3),3)
X <- num_ex["t=2"]
X[3,1] <- "--"
X[7,1] <- "--"
reg2 <- cbind(Y,X)
reg2 <- rename(reg2,"Y"="t=3","X"="t=2")

kable(reg2,caption="Regression at time 2")

```

We then regress $Y$ on a constant $X$ and $X^2$ which in this case results in the function $E[Y|X]=-13.35+27.98X-14.59X^2$. With this simple conditional expectation model it is possible to compare the value of immediate exercise at time 2, with the value from continuation. 
```{r echo=FALSE, cache=TRUE}
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

```{r echo=FALSE, cache=TRUE}
empty1 <-  rename(as.data.frame(matrix("--",8,1)),"t=1"="V1")
t3 <- rename(as.data.frame(matrix("0",8,1)),"t=3"="V1")
e2 <- opt_early$Exercise
cf2 <- cbind(empty1,rename(cbind(e2,t3),"t=2"="e2"))
cf2[3,2]=0
cf2[7,2]=0
kable(cf2, caption="Cash-flow matrix at time 2")
```

By recursively repeating the process for each time step we get the single optimal stopping time for each simulation path. It is then possible to average over all paths and compute the value of the American option. 
