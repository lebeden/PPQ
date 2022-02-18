### Average likelihood experiment
library(PPQ)
library(tidyverse)
options(max.print = 6)
# Step 0: set parameters --------------------------------------------------

eta <- 1
mu <- 1
s <- 20
# EB (Job size epxectation)eta/mu
n <- 1000
gamma <- 10
lambda_0 <- 10
theta <- 2.5e
PARAMS <-c(gamma=gamma,lambda_0=lambda_0,theta=theta)
(rho <- (gamma + lambda_0) * eta / (s*mu))


# Step 1: Generate many datasets ----------------------------------------------
MC <- 100
res_list <- list()
tic()
for(i in 1:MC){
  set.seed(i)
  RES <- resSimCosine(n=n,gamma = gamma,lambda_0 = lambda_0,theta = theta,s = s,eta = eta,mu = mu)

  A <- RES$A
  W <- RES$Wj
  X <- RES$Xj

  dati <- data.frame(A,W,X)
  res_list[[i]] <- dati
  svMisc::progress(i,MC)
}


# Step 2 - obtain likelihood function from each sample --------------------

makeLikelihoodFunction <- function(dati){

  # the i-th likelihood function
  ithLL <- function(gamma,lambda_0,theta){
    A <- dati$A
    W <- dati$W
    X <- dati$X
    A_i = A[-1]
    A_tilde <-c(0,cumsum(A))
    A_tilde <- A_tilde[-length(A_tilde)]
    A_tilde_i = cumsum(A_i)
    W_i = W[-1]
    w_i = W[-length(W)]
    x_i = X[-length(X)]

    logLikelihood <-
      log(gamma/2 + lambda_0 + (gamma*cos(2*pi*A_tilde_i))/2) +
      log(exp(-W_i*theta)) +
      (gamma*exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/
      (2*(4*pi^2 + theta^2)) - (lambda_0*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/theta - (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta)

    return(-mean(logLikelihood))
  }

  return(ithLL)
}

fun_list <- lapply(res_list,makeLikelihoodFunction)


pointAverageNegLogLikelihood <- function(gamma,lambda_0,theta,fun_list){
  params <- list(gamma,lambda_0,theta)
  sapply(fun_list, function(f) do.call(f,params) ) %>%
    mean
}


pointAverageNegLogLikelihood(fun_list = fun_list,gamma = 10, lambda_0 = 10, theta = 2)

gridAverageNegLogLikelihood <- Vectorize(pointAverageNegLogLikelihood,
                                         vectorize.args = c("gamma","lambda_0","theta"))

g <- 10:20
l <- 10:20
gr <- expand.grid(g,l)
th <- 2
a <- gridAverageNegLogLikelihood(fun_list=fun_list,gamma = gr$Var1,lambda_0=gr$Var2,th)

gamma.grid <- seq(gamma/2,gamma*3/2,length.out = 10)
lambda_0.grid <- seq(lambda_0/2,lambda_0*3/2,length.out = 10)
theta.grid <- seq(theta/2,theta*3/2,length.out = 10)
grids <- expand.grid(gamma.grid,lambda_0.grid,theta.grid) %>% as.data.frame()
names(grids) <- c("gamma","lambda_0","theta")
mapply(pointAverageNegLogLikelihood,
       fun_list = fun_list,
       gamma = grids$gamma,
       lambda_0 = grids$lambda_0,
       theta = grids$theta)
library(rgl)
open3d()
plot.function()
wire3d(mesh3d(g,l,a))

negLogLikelihoodMean()
# Functions ---------------------------------------------------------------

# this is an auxiliary function to compute the loglikelihood from the stored dataframes
# params = c(gamma,lambda_0, theta)
datLL <- function(dati, gamma, lambda_0, theta){
  A <- dati$A
  W <- dati$W
  X <- dati$X
  A_i = A[-1]
  A_tilde <-c(0,cumsum(A))
  A_tilde <- A_tilde[-length(A_tilde)]
  A_tilde_i = cumsum(A_i)
  W_i = W[-1]
  w_i = W[-length(W)]
  x_i = X[-length(X)]

  logLikelihood <-
    log(gamma/2 + lambda_0 + (gamma*cos(2*pi*A_tilde_i))/2) +
     log(exp(-W_i*theta)) +
    (gamma*exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/
     (2*(4*pi^2 + theta^2)) - (lambda_0*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/theta - (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta)

  return(-sum(logLikelihood))
}

