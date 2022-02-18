# testing file
library(tictoc)

library(simsalapar)
library(parallel)
library(PPQ)
N_SIM <- 80
# frozen parameters
LAMBDA_0 <- 10
ETA <- 1
MU <- 1
# grid parameters:
SAMPLE_SIZES <- 5e4
GAMMA <- c(3,10,19)
THETA <- c(1,2.5,5)
SERVERS <- c(5,10,20)

vList <- varlist(
  n.sim = list(value = N_SIM), # type = N
  n     = list(type="grid", value = SAMPLE_SIZES), # sample sizes
  lambda_0 = list(type = "frozen", value = LAMBDA_0),
  gamma = list(type = "grid",value = GAMMA),
  theta = list(type = "grid", value = THETA),
  s = list(type = "grid", value = SERVERS),
  eta = list(type = "frozen", value = ETA),
  mu = list(type = "frozen", value = MU)
)
cl <- makeCluster(7)
tic()
a1 <- doForeach(vList = vList, doOne = atomicSim,extraPkgs = c("PPQ"),cluster = cl)
toc()
## note the object in environment!
## a is res for nsim = 8
## a1 is for nsim = 80
# simulate ----------------------------------------------------------------

tic()
eta <- 1
mu <- 1
s <- 20
# EB (Job size epxectation)eta/mu
n <- 50000
gamma <- 10
lambda_0 <- 10
theta <- 2.5
(rho <- (gamma + lambda_0) * eta / (s*mu))

PARAMS <- c(gamma = gamma,lambda_0=lambda_0,theta = theta)
tic()
RES <- resSimCosine(n=n,gamma = gamma,lambda_0 = lambda_0,theta = theta,s = s,eta = eta,mu = mu)
toc()
# estimate ----------------------------------------------------------------
RES$Pl # % customers lost
rho # max load

PARAMS # true params
tic()
(boris <- mleBoris(RES,PARAMS,type = "s"))
(liron <- mleLironThetaLambda(RES))
toc()



opt <-
  optim(PARAMS, # note that PARAMS is temporary
        fn = negLogLikelihoodFull,
        lower = PARAMS/2,
        upper = PARAMS*2,
        method = "L-BFGS-B",
        gr = gradNegLogLikelihoodFull,
        RES = RES)
opt
opt$par
RES$Pl
negLogLikelihoodFull(c(10,10,1),RES)
gradNegLogLikelihoodFull(opt$par,RES)

grid_opt <- NMOF::gridSearch(negLogLikelihoodFull,
                 lower = PARAMS/2,
                 upper = PARAMS*2,
                 n = 20,

                 RES = RES
                 )
grid_opt$minlevels

dfoptim::nmkb(PARAMS, # note that PARAMS is temporary
              fn = negLogLikelihoodFull,
              lower = PARAMS/2,
              upper = PARAMS*2,RES=RES)

minqa::bobyqa(PARAMS, # note that PARAMS is temporary
              fn = negLogLikelihoodFull,
              lower = PARAMS/2,
              upper = PARAMS*2,RES=RES,)

grid_opt$minlevels



# -------------------------------------------------------------------------


