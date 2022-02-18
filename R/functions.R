# Simulation --------------------------------------------------------------


#' Rate function factory
#'
#' @param gamma
#' @param lambda_0
#'
#' @return A function with argument t that computes the arrival rate
#' @export
#'
#' @examples
#' rate <- RATE(5,10)
#' rate(1:3)
RATE <- function(gamma,lambda_0){
  return(function(t) lambda_0 + (gamma/2) * (1 + cos(2*pi*t)))
}




#' Next Arrival of inhomogeneous Poisson process
#'
#' Generates the next arrival time (not inter-arrival!) of a nonhomogeneous Poisson process.
#' @param current_time The current time (a non-negative real number)
#' @param gamma Parameter of the cosine arrival
#' @param lambda_0 Parameter of the cosine arrival
#' @return The time of the next arrival.
#'
#' @export
nextArrivalCosine <- function(current_time, gamma, lambda_0){
  lambda_sup <- gamma + lambda_0 # highest rate in the cosine model
  rate <- RATE(gamma = gamma, lambda_0 = lambda_0)
  arrived <- FALSE
  while (!arrived) {
    u1 <- runif(1)
    current_time <- current_time - (1/lambda_sup) * log(u1) # generate next arrival from Pois(sup_lambda)
    u2 <- runif(1)
    arrived <- u2 <= rate(current_time) / lambda_sup
  }
  return(current_time)
}



#' Customer's job size and patience
#' Provides with a realization of Y~Exp(theta) and B~Gamma(eta,mu)
#' @param theta The parameter of the exponential patience.
#' @param eta Shapae parameter of the job size.
#' @param mu rate parameter of the job size.
#'
#' @return A list with elements 'Patience' and 'Jobsize'
#' @export
#'
#' @examples
customerExp <- function(theta,eta,mu){

  B <-  rgamma(1, shape = eta, rate = mu) #Job sizes
  Y <- rexp(1,theta)
  return(list(Patience=Y,Jobsize=B))
}



#' Liron's Virtual Waiting function
#'
#' @param Res.service the vector of residual service times
#' @param s the number of servers
#' @export
VW <- function(Res.service,s){
  Q.length <- length(Res.service) #Number in the system
  if(Q.length < s){virtual.wait <- 0}else
  {
    D.times <- rep(NA,Q.length+1)
    Vw <- rep(NA,Q.length+1)
    D.times[1:s] <- Res.service[1:s]
    Vw[1:s] <- rep(0,s) #VW for customers in service
    for(i in (s+1):(Q.length+1))
    {
      D.i <- sort(D.times[1:i]) #Sorted departures of customers ahead of i
      Vw[i] <- D.i[i-s] #Virtual waiting time for position i
      if(i <= Q.length){D.times[i] <- Res.service[i]+Vw[i]} #Departure times
    }
    virtual.wait <- Vw[Q.length+1]
  }
  return(virtual.wait)
}

#' Atomic Simulation for periodic arrivals (cosine model)
#' @param gamma The periodic arrival rate maximum amplitude
#' @param lambda_0 The constant arrival rate
#' @param theta Parameter of the exponential patience
#' @param eta Shape parameter of job size.
#' @param mu Rate parameter of job size.
#' @param s Number of servers.
#' @param n Number of samples to generate.
#' @return A list with the estimates
#' @export
#' @examples
resSimCosine <- function(n, gamma, lambda_0, theta, s, eta, mu){
  #find the supremum of the arrival function
  lambda_sup <- lambda_0 + gamma
  lambdaFunction <-   RATE(gamma = gamma, lambda_0 = lambda_0)
  #Running simulation variables
  klok <- 0
  n.counter <- 0 #Observation counter
  Q.length <- 0 #Queue length (including service)
  virtual.wait <- 0 #Virtual waiting time process
  m <- 0 #Total customer counter
  #A <- rexp(1,lambda) #First arrival time
  A <- nextArrivalCosine(klok[length(klok)],  gamma, lambda_0)
  A <- A - klok[length(klok)]
  klok <- c(klok, klok[length(klok)] + A)
  trans.last <- A #Counter of time since last transition
  time.last <- A #Counter of time since last admission event
  Res.service <- numeric(0) #Vector of residual service times

  #Output vectors:
  Wj <- rep(NA,n) #Vector of workloads before jumps
  Xj <- rep(NA,n) #Vector of workload jumps
  Aj <- rep(NA,n) #Vector of effective inter-arrival times
  Qj <- rep(NA,n) #Vector of queue lengths at arrival times
  IPj <- A #Vector of idle periods
  Yj <- rep(NA,n) # Vector of patience values of customers that join
  Q.trans <- 0 #Vector of queue lengths at transitions
  IT.times <- numeric(0) #Vector of inter-transition times
  Pl <- 1 #Proportion of lost customers
  Nl <- 0 # number of lost customers

  while(n.counter < n+1)
  {
    m <- m+1
    customer <- customerExp(theta = theta, eta = eta, mu = mu)
    #Generate patience and job size:
    B <- customer$Jobsize
    Y <- customer$Patience
    if(virtual.wait <= Y) #New customer if patience is high enough
    {

      n.counter <- n.counter+1 #Count observation
      Res.service <- c(Res.service,B) #Add job to residual service times vector
      Q.length <- length(Res.service) #Queue length
      Q.trans <- c(Q.trans,Q.length) #Add current queue length
      #Add new observations:
      Wj[n.counter] <- virtual.wait
      Aj[n.counter] <- time.last
      Qj[n.counter] <- Q.length-1 #Queue length (excluding new arrival)
      Yj[n.counter] <- Y # patience of the customer arriving
      IT.times <- c(IT.times,trans.last) #Update transition time
      trans.last <- 0 #Reset transition time
      time.last <- 0 #Reset last arrival time
      Pl <- m*Pl/(m+1) #Update loss proportion
    }else { Pl <- m*Pl/(m+1)+1/(m+1)
    Nl <- Nl + 1}

    #Update system until next arrival event
    #A <- rexp(1,lambda) #Next arrival time
    A <- nextArrivalCosine(klok[length(klok)],  gamma, lambda_0)
    A <- A - klok[length(klok)]
    klok <- c(klok, klok[length(klok)] + A)
    time.last <- time.last+A #Add arrival time to effective arrival time

    #Departure and residual service times of customers in the system:
    Q.length <- length(Res.service) #Queue length
    D.times <- rep(NA,Q.length) #Departure times
    Vw <- rep(NA,Q.length) #Virtual waiting times
    for(i in 1:Q.length)
    {
      if(i <= s)
      {
        Vw[i] <- 0 #No virtual waiting time
        D.times[i] <- Res.service[i] #Departure time is the residual service time
        Res.service[i] <- max(Res.service[i]-A,0) #Update residual service time
      }else
      {
        D.i <- sort(D.times[1:i]) #Sorted departures of customers ahead of i
        Vw[i] <- D.i[i-s] #Time of service start for customer i
        D.times[i] <- Res.service[i]+Vw[i] #Departure time
        serv.i <- max(0,A-Vw[i]) #Service obtained before next arrival
        Res.service[i] <- max(Res.service[i]-serv.i,0) #New residual service
      }
    }
    #Jump of virtual waiting time:
    if(virtual.wait <= Y)
    {
      if(Q.length < s)
      {
        Xj[n.counter] <- 0
      }else
      {
        Xj[n.counter] <- sort(D.times)[Q.length+1-s]-virtual.wait
      }
    }
    #Update residual service times:
    Res.service <- Res.service[!(Res.service == 0)] #Remove completed services

    #Update transition times and queue lengths:
    D.before <- which(D.times <= A) #Departing customers before next arrival
    if(length(D.before) > 0)
    {
      T.d <- sort(D.times[D.before]) #Sorted departure times
      for(i in 1:length(D.before))
      {
        Q.trans <- c(Q.trans,Q.length-i) #Update queue length at departures
        if(i ==1)
        {
          trans.last <- trans.last+T.d[1] #Update time since last transition
          IT.times <- c(IT.times,trans.last) #Departure transition time
          trans.last <- 0 #Reset transition time
        }else
        {
          trans.last <- trans.last+T.d[i]-T.d[i-1] #Update time since last transition
          IT.times <- c(IT.times,trans.last) #Departure transition time
          trans.last <- 0 #Reset transition time
        }
        if(Q.trans[length(Q.trans)] == 0){IPj <- cbind(IPj,A-T.d[length(T.d)])} #Add idle time observation
      }
      trans.last <- A-T.d[i] #Update remaining time until next arrival
    }else if (length(D.before) == 0 )
    {
      trans.last <- trans.last+A #Update timer since least transition with new arrival
    }
    virtual.wait <- VW(Res.service,s) #Update virtual waiting time

  }
  # progress bar
  if (m%%1000 == 0 )
    cat("* ")
  RES <- list(Aj=Aj,Xj=Xj,Wj=Wj,Qj=Qj,
              IPj=IPj,Q.trans=Q.trans,IT.times=IT.times,
              Yj = Yj, klok = klok, Pl = Pl, Nl = Nl)

  return(RES)
}



#' Function to make many files with simulation results
#'
#' @param N_files Desired number of files to be created in the working directory
#' @param n_obs Number of effective arrivals per file. Defaults to 10,000. Do not use less than 5000.
#' @param gamma The periodic arrival rate maximum amplitude
#' @param lambda_0 The constant arrival rate
#' @param theta Parameter of the exponential patience
#' @param eta Shape parameter of job size.
#' @param mu Rate parameter of job size.
#' @param s Number of servers.
#' @details The filenames have the parameter values encoded alongside a timestamp,
#' which is meant for protection against overwriting in case of repeated calls to this function.
#' Also note that if files take less than a second to generate, they will not be written.
#' @return
#' @export
#'
#' @examples
makeSimFilesAWX <- function(N_files, n_obs, gamma, lambda_0, theta, s, eta, mu){
  for (i in 1:N_files){
    RES <- resSimCosine(n=n_obs,gamma = gamma,lambda_0 = lambda_0,theta = theta,s = s,eta = eta,mu = mu)
    A <- RES$A
    W <- RES$Wj
    X <- RES$Xj
    dat <- data.frame(A=A,W=W,X=X)
    name <- filenamer::filename(paste0("AWX_","n=",n_obs),
                                ext = "csv",subdir = FALSE)
    param_values <- c("gamma"=gamma,"lambda_0"=lambda_0,"theta"=theta,"s"=s,"eta"=eta,"mu" = mu)
    # making the tag - use only first letter of each parameter name
    totag <- paste0(substr(names(param_values),1,1),"=",param_values)
    name <- insert(name,totag)
    name <- as.character(name)
    write.csv(dat,file = name,row.names = FALSE)

  }
}

# Estimation --------------------------------------------------------------


#' (negative) log-likelihood for the sinusoidal arrivals model
#'
#' @param params A vector with values (gamma,lambda_0, theta).
#' @param RES simulation results
#' @return The negative log-likelihood at the point provided.
#' @export
#'
negLogLikelihoodFull <- function(params,RES){
  gamma <- params[1]
  lambda_0 <- params[2]
  theta <- params[3]
  # params is a vector of (lambda_bar,alpha,theta)
  W <- RES$Wj
  Q.trans <- RES$Q.trans
  IT.times <- RES$IT.times
  Transitions <- c(0,cumsum(IT.times))
  A <- RES$A
  A_tilde <-c(0,cumsum(A))
  A_tilde <- A_tilde[-length(A_tilde)]
  X <- RES$Xj
  Q <- RES$Qj
  Pl <- RES$Pl
  Y <- RES$Yj
  # data from RES
  A_i = A[-1];
  A_tilde_i = cumsum(A_i);
  W_i = W[-1]
  w_i = W[-length(W)];
  x_i = X[-length(X)]

  # obtained by symbolic calculation in MATLAB
  logLikelihood <-
  log(gamma/2 + lambda_0 + (gamma*cos(2*pi*A_tilde_i))/2) +
    log(exp(-W_i*theta)) +
    (gamma*exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/
    (2*(4*pi^2 + theta^2)) - (lambda_0*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/theta - (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta)

  negLL <- - sum(logLikelihood)

  return(negLL) # NEGATIVE LOG-LIKELIHOOD
}

#' (negative) mean log-likelihood for the sinusoidal arrivals model
#' The mean is returned instead of the sum - should help gradient based optimizers
#' @param params A vector with values (gamma,lambda_0, theta).
#' @param RES simulation results(should be NULL if dati supplied)
#' @param dati compact simulation results for memory economy (A,W,X).
#' @return The negative log-likelihood at the point provided.
#' @export
#'
negLogLikelihoodMean<- function(params,RES,dati){

  # params is a vector of (lambda_bar,alpha,theta)
  if(!is.null(dati)){
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


  }
  else if (is.null(dati) && !is.null(RES)){
    W <- RES$Wj
    Q.trans <- RES$Q.trans
    IT.times <- RES$IT.times
    Transitions <- c(0,cumsum(IT.times))
    A <- RES$A
    A_tilde <-c(0,cumsum(A))
    A_tilde <- A_tilde[-length(A_tilde)]
    X <- RES$Xj
    Q <- RES$Qj
    Pl <- RES$Pl
    Y <- RES$Yj
    # data from RES
    A_i = A[-1];
    A_tilde_i = cumsum(A_i);
    W_i = W[-1]
    w_i = W[-length(W)];
    x_i = X[-length(X)]

  }
  else stop("no data supplied!")


  # obtained by symbolic calculation in MATLAB
  logLikelihood <-
    log(gamma/2 + lambda_0 + (gamma*cos(2*pi*A_tilde_i))/2) +
    log(exp(-W_i*theta)) +
    (gamma*exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/
    (2*(4*pi^2 + theta^2)) - (lambda_0*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/theta - (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta)

  # mean instead of sum, the minus is for the negative log-likelihood
  negLL <- - mean(logLikelihood)

  return(negLL) # NEGATIVE LOG-LIKELIHOOD
}

#' Wrapper for negative log likelihood which accepts separate parameters
#'
#' @param gamma
#' @param lambda_0
#' @param theta
#' @param RES
#' @param dati
#'
#' @return
#' @export
#'
#' @examples
nll.mean <- function(gamma,lambda_0,theta,RES,dati){
  params <- numeric(3)
  params[1] <-   gamma
  params[2] <- lambda_0
  params[3] <- theta
  if(!is.null(dati)){
    negLogLikelihoodMean(params = params,dati = dati)
  }
  else if (!is.null(RES) )

}

gradNegLogLikelihoodFull <- function(params,RES){
  gamma <- params[1]
  lambda_0 <- params[2]
  theta <- params[3]
  # params is a vector of (lambda_bar,alpha,theta)
  W <- RES$Wj
  Q.trans <- RES$Q.trans
  IT.times <- RES$IT.times
  Transitions <- c(0,cumsum(IT.times))
  A <- RES$A
  A_tilde <-c(0,cumsum(A))
  A_tilde <- A_tilde[-length(A_tilde)]
  X <- RES$Xj
  Q <- RES$Qj
  Pl <- RES$Pl
  Y <- RES$Yj
  # data from RES
  A_i = A[-1];
  A_tilde_i = cumsum(A_i);
  W_i = W[-1]
  w_i = W[-length(W)];
  x_i = X[-length(X)]

  # obtained by symbolic calculation in MATLAB
  # derivative wrt gamma
  dg <- (cos(2*pi*A_tilde_i)/2 + 1/2)/(gamma/2 + lambda_0 + (gamma*cos(2*pi*A_tilde_i))/2) - (exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta) + (exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(2*(4*pi^2 + theta^2))
  # derivative wrt lambda_0
  dl <- log(gamma/2 + lambda_0 + (gamma*cos(2*pi*A_tilde_i))/2) + log(exp(-W_i*theta)) + (gamma*exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(2*(4*pi^2 + theta^2)) - (lambda_0*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/theta - (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta)
  # derivative wrt theta
  dt <- (lambda_0*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/theta^2 - W_i - (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i)) - cos(2*pi*A_tilde_i) + 2*A_i*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) + A_i*theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(2*(4*pi^2 + theta^2)) + (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta^2) + (gamma*exp(-theta*(w_i + x_i))*(w_i + x_i)*(exp(A_i*theta) - 1))/(2*theta) - (gamma*exp(-theta*(w_i + x_i))*(w_i + x_i)*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(2*(4*pi^2 + theta^2)) + (lambda_0*exp(-theta*(w_i + x_i))*(w_i + x_i)*(exp(A_i*theta) - 1))/theta - (gamma*theta*exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(4*pi^2 + theta^2)^2 - (A_i*gamma*exp(A_i*theta)*exp(-theta*(w_i + x_i)))/(2*theta) - (A_i*lambda_0*exp(A_i*theta)*exp(-theta*(w_i + x_i)))/theta
  # the NEGATIVE gradient vector (hence the minus):
  g <- -1 * c(sum(dg), sum(dl), sum(dt))

  return(g)


}

gradNegLogLikelihoodMean <- function(params,RES){
  gamma <- params[1]
  lambda_0 <- params[2]
  theta <- params[3]
  # params is a vector of (lambda_bar,alpha,theta)
  W <- RES$Wj
  Q.trans <- RES$Q.trans
  IT.times <- RES$IT.times
  Transitions <- c(0,cumsum(IT.times))
  A <- RES$A
  A_tilde <-c(0,cumsum(A))
  A_tilde <- A_tilde[-length(A_tilde)]
  X <- RES$Xj
  Q <- RES$Qj
  Pl <- RES$Pl
  Y <- RES$Yj
  # data from RES
  A_i = A[-1];
  A_tilde_i = cumsum(A_i);
  W_i = W[-1]
  w_i = W[-length(W)];
  x_i = X[-length(X)]

  # obtained by symbolic calculation in MATLAB
  # derivative wrt gamma
  dg <- (cos(2*pi*A_tilde_i)/2 + 1/2)/(gamma/2 + lambda_0 + (gamma*cos(2*pi*A_tilde_i))/2) - (exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta) + (exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(2*(4*pi^2 + theta^2))
  # derivative wrt lambda_0
  dl <- log(gamma/2 + lambda_0 + (gamma*cos(2*pi*A_tilde_i))/2) + log(exp(-W_i*theta)) + (gamma*exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(2*(4*pi^2 + theta^2)) - (lambda_0*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/theta - (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta)
  # derivative wrt theta
  dt <- (lambda_0*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/theta^2 - W_i - (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i)) - cos(2*pi*A_tilde_i) + 2*A_i*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) + A_i*theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(2*(4*pi^2 + theta^2)) + (gamma*exp(-theta*(w_i + x_i))*(exp(A_i*theta) - 1))/(2*theta^2) + (gamma*exp(-theta*(w_i + x_i))*(w_i + x_i)*(exp(A_i*theta) - 1))/(2*theta) - (gamma*exp(-theta*(w_i + x_i))*(w_i + x_i)*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(2*(4*pi^2 + theta^2)) + (lambda_0*exp(-theta*(w_i + x_i))*(w_i + x_i)*(exp(A_i*theta) - 1))/theta - (gamma*theta*exp(-theta*(w_i + x_i))*(2*pi*sin(2*pi*A_tilde_i) + theta*cos(2*pi*A_tilde_i) - 2*pi*sin(2*pi*(A_i + A_tilde_i))*exp(A_i*theta) - theta*exp(A_i*theta)*cos(2*pi*(A_i + A_tilde_i))))/(4*pi^2 + theta^2)^2 - (A_i*gamma*exp(A_i*theta)*exp(-theta*(w_i + x_i)))/(2*theta) - (A_i*lambda_0*exp(A_i*theta)*exp(-theta*(w_i + x_i)))/theta
  # the NEGATIVE gradient vector (hence the minus):
  g <- -1 * c(mean(dg), mean(dl), mean(dt))

  return(g)


}


#' MLE for the cosine arrival + exponential patience model
#'
#' @param RES
#' @param PARAMS (temporary) True values of parameters to use a starting points
#' @param type type of log-likelihood to be optimized - sum ("s") or mean ("m")
#' @return gradient vector of the negative log-likelihood
#' @export
#'
#' @examples
mleBoris <- function(RES, PARAMS ,type){
  if (type == "s"){
    opt <-
    optim(PARAMS, # note that PARAMS is temporary
          fn = negLogLikelihoodFull,
          lower = PARAMS/2,
          upper = PARAMS*2,
          method = "L-BFGS-B",
          gr = gradNegLogLikelihoodFull,
          RES = RES)
  }
  if (type == "m"){
    opt <-
      optim(PARAMS, # note that PARAMS is temporary
          fn = negLogLikelihoodMean,
          lower = PARAMS/2,
         upper = PARAMS*2,
         method = "L-BFGS-B",
          gr = gradNegLogLikelihoodMean,
          RES = RES)
  }
  ans <- opt$par
  names(ans) <- c("gamma","lambda_0","theta")
  return(ans)
}

# Liron's estimator -------------------------------------------------------


#Lambda MLE given an estimator for theta (exponential patience)
lambda.MLE <- function(theta.hat,A,W,X){
  n <- length(W)
  a <- exp(-theta.hat*W[2:n])-exp(-theta.hat*(W[1:(n-1)]+X[1:(n-1)]))
  b <- theta.hat*pmax(rep(0,n-1),A[2:n]-W[1:(n-1)]-X[1:(n-1)])
  lambda.hat <- n*theta.hat/sum(a+b)
  return(lambda.hat)
}

mleLironThetaLambda<- function(RES,acc=1e-4){
  A <- RES$Aj
  W <- RES$Wj
  X <- RES$Xj
  n <- length(W)
  Theta <- c(0,10^3) #Search range
  d <- acc*2
  #Bisection search for optimal theta
  while(abs(d)>acc)
  {
    theta.hat <- mean(Theta)
    lambda.hat <- lambda.MLE(theta.hat,A,W,X) #Lambda mle for given theta
    a <- (1+theta.hat*W[2:n])*exp(-theta.hat*W[2:n])-(1+theta.hat*(W[1:(n-1)]+X[1:(n-1)]))*exp(-theta.hat*(W[1:(n-1)]+X[1:(n-1)]))
    d <- mean(W[2:n])-lambda.hat*mean(a)/(theta.hat^2)
    #Update value:
    if(d > acc){Theta[2] <- theta.hat}
    if(d < -acc){Theta[1] <- theta.hat}
  }
  return(data.frame(theta.hat=theta.hat,lambda.hat=lambda.hat))
}


# Parallel Simulation -----------------------------------------------------

#' atomic realization for parallel simulation study
#'
#' Atomic Simulation for periodic arrivals (cosine model)
#' @param n Number of samples (effective arrivals) to generate.
#' @param lambda_0 The constant arrival rate
#' @param gamma The periodic arrival rate maximum amplitude
#' @param theta Parameter of the exponential patience
#' @param eta Shape parameter of job size.
#' @param mu Rate parameter of job size.
#' @param s Number of servers.

#'
#' @return A vector with Boris's etimates and Liron's etsimates
#' @export
#'
#' @examples
atomicSim <- function(n,lambda_0, gamma,theta, eta =1, mu = 1, s ){
  # temporary: use true parameters as starting values
  PARAMS <- c(gamma = gamma,lambda_0=lambda_0,theta = theta)
    # generate the sample
  RES <- resSimCosine(n=n,gamma = gamma,lambda_0 = lambda_0,theta = theta,s = s,eta = eta,mu = mu)
  boris <- mleBoris(RES,PARAMS)
  liron <- mleLironThetaLambda(RES)
  ans <- as.numeric(c(boris,liron))
  names(ans) <- c("B_gamma", "B_lambda", "B_theta", "L_theta", "L_lambda")
  return(ans)
}

