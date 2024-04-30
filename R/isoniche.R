
#' isoniche
#'
#' Function to fit a bivariate regression model for calculating standard ellipses
#' in isotopic niche comparisons.
#'
#' @param mean A list of length 2 with formulas defining the models for the means along
#' each axis.
#' @param var A list of length 3 of one-sided formulas that define the models for the
#' scale parameters (standard deviations) along each axis as well as a model for the correlation
#' parameter. These must be defined, but can be set to \code{~ 1} for constant scales or correlation
#' across groups / individuals.
#' @param data Data frame with one row per individual and columns for the isotope concentrations as
#' well as group designations (as factors) and any control variables.
#' @param ... Additional arguments that can be passed to the sampler. These include \code{iter}, \code{chains},
#' \code{thin}, \code{cores}, \code{control}. See \code{?rstan::sampling()} for more details.
#'
#' @return A fitted model object from the \pkg{rstan} package.
#'
#' @example modelfit_eg.R
#'
isoniche <- function(mean, var, data, ...){

  stan_vars <- list(...)
  if(is.null(stan_vars$iter)){
    stan_vars$iter <- 2000
  }
  if(is.null(stan_vars$chains)){
    stan_vars$chains <- 4
  }
  if(is.null(stan_vars$thin)){
    stan_vars$thin <- 1
  }
  if(is.null(stan_vars$cores)){
    stan_vars$cores <- 1
  }
  if(length(mean) != 2){
    stop("Incorrect number of model formulas specifying the mean. Expecting 2.\n")
  }
  if(length(var) != 3){
    stop("Incorrect number of model formulas specifying the variance.\n
         Expecting 3; one for each scale parameter and one for the correlation parameter.\n")
  }
  # first create model matrices
  Xs <- lapply(mean, model.matrix, data = data)
  Zs <- lapply(var, model.matrix, data = data)

  # bind these together for joint model specification
  X <- do.call(cbind, Xs)
  Z <- do.call(cbind, Zs[1:2])

  # final model matrix defines the model for the correlation
  G <- Zs[[3]]

  # extract y variables
  ynames <- sapply(mean, all.vars)[1, ]

  # compile variables for stan
  datlist <- list(
    N = nrow(X),
    P = sapply(Xs, ncol),
    K = sapply(Zs[1:2], ncol),
    J = ncol(G),
    X = X, Z = Z, G = G,
    y = as.matrix(data[, ynames])
  )

  rstan::sampling(
    stanmodels$isonich,
    data = datlist,
    iter = stan_vars$iter,
    chains = stan_vars$chains,
    thin = stan_vars$thin,
    cores = stan_vars$cores,
    control = stan_vars$control
  )

}

