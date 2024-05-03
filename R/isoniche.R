
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
#' @example /inst/examples/isoniche_eg.R
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

  fit <- rstan::sampling(
    stanmodels$isonich,
    data = datlist,
    iter = stan_vars$iter,
    chains = stan_vars$chains,
    thin = stan_vars$thin,
    cores = stan_vars$cores,
    control = stan_vars$control
  )

  return(
    list(
      fit = fit,
      model = list(
        mean = mean,
        var = var
      ),
      data = list(
        df = data,
        datlist = datlist
      )
    )
  )

}


#' Construct standard ellipses from fitted model
#'
#' This function takes in a dataframe that represents the new data over which
#' to predict. For example, a set of groups for which to construct standard
#' ellipses.
#'
#' @param mfit Fitted isoniche model object.
#' @param newdat New data, with all the same columns as the original data used to
#' fit the model.
#' @param n Number of draws from the posterior that can be used to display standard ellipses.
#' The default is 1, which results in ellipses being drawn from the marginal means of the posterior
#' of each parameter. If \code{n} > 1, then the function will draw \code{n} samples from the joint
#' posterior and construct a standard ellipse for each set of parameters.
#'
#' @return A dataframe that can be used for plotting standard ellipses or calculating isotopic niche
#' statistics.
#'
#' @examples /inst/examples/isoniche_eg.R
construct_ellipses <- function(mfit, newdat, n = 1){

  varnames <- as.vector(sapply(mfit$model$mean, all.vars)[-1, ])
  varnames <- c(
    varnames,
    as.vector(sapply(mfit$model$var, all.vars))
  )
  if(any(!(unique(varnames) %in% names(newdat)))){
    stop("newdat must contain columns for all the variables used to fit the model.\n")
  }
  # construct new model matrices for prediction
  Xs <- lapply(mfit$model$mean, model.matrix, data = newdat)
  Zs <- lapply(mfit$model$var, model.matrix, data = newdat)

  # bind these together for joint model specification
  X_new <- do.call(cbind, Xs)
  Z_new <- do.call(cbind, Zs[1:2])

  # final model matrix defines the model for the correlation
  G_new <- Zs[[3]]

  # construct parameter matrices
  if(n == 1){

  }

}


# helper functions ---------------------------


make_parmat <- function(P, theta_1, theta_2){

  B <- matrix(data = 0, nrow = sum(P), ncol = 2)
  B[1:P[1], 1] <- theta_1
  B[(P[1] + 1):sum(P), 2] <- theta_2

}





