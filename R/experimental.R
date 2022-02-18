
#' Make Parameter Grid for testing
#'
#' @param params named list of the parameters (gamma,lambda_0,theta)
#' @param spans vector of length 3 which determines the span of each axis in the
#' grid. Values must be in the interval [0,1). See details.
#' @param grid.size an integer vector of length 1 or 3. In the former case, all
#' axes have the same number of points, and in the latter each parameter is
#' matched to the corresponding element of `grid.size`.
#'
#' @return A dataframe with the parameter grids in the columns.
#' @export
#'
#' @examples
makeParGrid <- function(params,spans,grid.size){

  if (!length(spans) %in% c(1,3))
    stop("Wrong length of spans. Must be of lengths either 1 or 3.")
  if (any(spans >= 1 | spans <=0))
    stop("All parameters must be positive. Please correct the value of spans.")
  if (length(params) != 3)
    stop("Invalid number of parameters, must be three - gamma, lambda_0, theta.")
  # the central parameter values
  pars <- with(params,c(gamma, lambda_0, theta))
  grid.values <- mapply(seq,
                        from = pars * (1 - spans),
                        to = pars * (1 + spans),
                        length.out = grid.size) %>%
    as.data.frame()


  names(grid.values) <- c("gamma", "lambda_0", "theta")
  outer(grid.values$gamma,grid.values$lambda_0,l1,theta = 1)

  grid <- expand.grid(grid.values$gamma,grid.values$lambda_0,grid.values$theta)
  names(grid) <- c("gamma", "lambda_0", "theta")
  return(grid)
}
names(params) <- c("gamma", "lambda_0", "theta")
grid <- makeParGrid(params = params,spans = .4,grid.size = 10)

sap


outerGrid <- function(logL_i,
                      params,
                      which.constant,
                      spans = c(0.5,0.5),
                      grid.size = 10){

  constant.index <- which(names(params) == which.constant)

  if (which.constant == "gamma"){
    pars <- params[-constant.index]
    grid.values <- mapply(seq,
                          from = pars * (1 - spans),
                          to = pars * (1 + spans),
                          length.out = grid.size) %>%
      as.data.frame()
    names(grid.values) <- names(params)[-constant.index]
    outer(grid.values[,1],grid.values[,2],logL_i,gamma = params[["gamma"]])
  }
  pars <- params[-constant.index]
  grid.values <- mapply(seq,
                        from = pars * (1 - spans),
                        to = pars * (1 + spans),
                        length.out = grid.size) %>%
    as.data.frame()
  gr <- expand.grid(grid.values)

  outer(gr[,1],gr[,2],FUN = logL_i, )
  return(eval(parse(text = which.constant)))
}
undebug(outerGrid)
outerGrid(params = PARAMS,which.constant = "theta")



library(MASS)
# from the fitdistr example
set.seed(123)
x <- rgamma(100, shape = 5, rate = 0.1)
fit <- fitdistr(x, dgamma, list(shape = 1, rate = 0.1), lower = 0.001)
loglik <- function(shape, rate) sum(dgamma(x, shape=shape, rate=rate,
                                           log=TRUE))
loglik <- Vectorize(loglik)
xlim <- fit$estimate[1]+4*fit$sd[1]*c(-1,1)
ylim <- fit$estimate[2]+4*fit$sd[2]*c(-1,1)

mfrow3d(1, 2, sharedMouse = TRUE)
persp3d(loglik,
        xlim = xlim, ylim = ylim,
        n = 30)
zlim <- fit$loglik + c(-qchisq(0.99, 2)/2, 0)
next3d()
persp3d(loglik,
        xlim = xlim, ylim = ylim, zlim = zlim,
        n = 30)
