
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
#' @export
#' @examples /inst/examples/isoniche_eg.R
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
#' @importFrom stats formula model.matrix
#' @importFrom mvtnorm rmvnorm
#' @export
#' @examples /inst/examples/isoniche_eg.R
#'
construct_ellipses <- function(mfit, newdat, n = 1){

  resp_names <- as.vector(sapply(mfit$model$mean, all.vars)[1, ])
  varnames <- as.vector(sapply(mfit$model$mean, all.vars)[-1, ])
  varnames <- c(
    varnames,
    as.vector(sapply(mfit$model$var, all.vars))
  )
  if(is.list(varnames)){
    tokeep <- !sapply(varnames, function(x, a = character(0)){identical(x, a)})
    varnames <- varnames[tokeep]
  }
  if(any(!(unique(varnames) %in% names(newdat)))){
    stop("newdat must contain columns for all the variables used to fit the model.\n")
  }
  # construct new model matrices for prediction
  rhs_mforms <- lapply(mfit$model$mean, rhs_form)
  Xs <- lapply(rhs_mforms, model.matrix, data = newdat)
  Zs <- lapply(mfit$model$var, model.matrix, data = newdat)

  # bind these together for joint model specification
  X_new <- do.call(cbind, Xs)
  Z_new <- do.call(cbind, Zs[1:2])

  # final model matrix defines the model for the correlation
  G_new <- Zs[[3]]

  # construct matrix of xy coords for a unit circle
  xy_circ <- rbind(
    cos(seq(0, 2 * pi, length.out = 100)),
    sin(seq(0, 2 * pi, length.out = 100))
  )
  ones <- matrix(data = 1, 2, 100)

  if(n == 1){
    # construct parameter matrices
    B <- make_parmat(
      mfit$data$datlist$P,
      rstan::summary(mfit$fit, pars = "beta_1")$summary[, "mean"],
      rstan::summary(mfit$fit, pars = "beta_2")$summary[, "mean"]
    )
    Zeta <- make_parmat(
      mfit$data$datlist$K,
      rstan::summary(mfit$fit, pars = "zeta_1")$summary[, "mean"],
      rstan::summary(mfit$fit, pars = "zeta_2")$summary[, "mean"]
    )
    gamma <- rstan::summary(mfit$fit, pars = "gamma")$summary[, "mean"]

    # make means and lower cholesky factors for each row of newdat
    mu_new <- lapply(
      1:nrow(newdat),
      function(i, X_new, B){
        as.double(X_new[i, ] %*% B)
      },
      X_new, B
    )
    L_new <- lapply(
      1:nrow(newdat),
      function(i, Z_new, Zeta, G_new, gamma){
        s <- exp(as.double(Z_new[i, ] %*% Zeta))
        rho <- 2 * stats::plogis(as.double(G_new[i, ] %*% gamma)) - 1
        return(make_L2d(rho, s))
      },
      Z_new, Zeta, G_new, gamma
    )

    # make dfs for ellipses coords
    dfs_new <- ellipse_df_list(mu_new, L_new, ones, xy_circ)

    # make sure the names line up
    dfs_new <- lapply(dfs_new, `colnames<-`, resp_names)

    # construct one big dataframe
    df_full <- as.data.frame(newdat[rep(1:nrow(newdat), each = 100), ])
    names(df_full) <- names(newdat)
    df_full$ellipse_id <- factor(rep(1:nrow(newdat), each = 100))
    return(cbind(
      df_full,
      Reduce(rbind, dfs_new)
    ))
  }
  # now if there are multiple draws from posterior
  if(n > 1){
    # construct lists of parameter matrices
    post_draws <- rstan::extract(mfit$fit)
    draws <- sample(1:length(post_draws$lp__), size = n)
    B_list <- lapply(
      draws,
      function(i, P, post_draws){
        beta_1 <- as.double(post_draws$beta_1[i, ])
        beta_2 <- as.double(post_draws$beta_2[i, ])
        make_parmat(P, beta_1, beta_2)
      },
      P = mfit$data$datlist$P,
      post_draws = post_draws
    )
    Zeta_list <- lapply(
      draws,
      function(i, P, post_draws){
        zeta_1 <- as.double(post_draws$zeta_1[i, ])
        zeta_2 <- as.double(post_draws$zeta_2[i, ])
        make_parmat(P, zeta_1, zeta_2)
      },
      P = mfit$data$datlist$K,
      post_draws = post_draws
    )
    gamma_list <- lapply(
      draws,
      function(i, post_draws){as.double(post_draws$gamma[i, ])},
      post_draws
    )

    # for each row of newdat and each set of params, create a mean and
    # lower Cholesky factor
    mu_new <- lapply(
      B_list,
      function(B, X_new){
        lapply(
          1:nrow(X_new),
          function(i, X_new, B){
            as.double(X_new[i, ] %*% B)
          },
          X_new = X_new, B = B
        )
      },
      X_new = X_new
    )
    L_new <- mapply(
      function(Z, Zeta, G, gamma, newdat){
        lapply(
          1:nrow(newdat),
          function(i, Z, Zeta, G, gamma){
            s <- exp(as.double(Z[i, ] %*% Zeta))
            rho <- 2 * stats::plogis(as.double(G[i, ] %*% gamma)) - 1
            return(make_L2d(rho, s))
          },
          Z, Zeta, G, gamma
        )
      },
      Zeta_list, gamma_list,
      MoreArgs = list(
        Z = Z_new, G = G_new, newdat = newdat
      ),
      SIMPLIFY = F
    )

    # Now, for each set of params, create a list of dfs
    # defining the plotting dataframes for each row of newdat
    dfs_list <- mapply(
      ellipse_df_list,
      mu_new, L_new,
      MoreArgs = list(
        ones, xy_circ
      ),
      SIMPLIFY = F
    )

    # make sure column names align
    dfs_list <- lapply(
      dfs_list,
      function(x, resp_names){
        lapply(x, `colnames<-`, resp_names)
      },
      resp_names
    )

    # construct one giant df
    row_reps <- rep(rep(1:nrow(newdat), each = 100), n)
    df_full <- as.data.frame(newdat[row_reps, ])
    names(df_full) <- names(newdat)
    df_full$draw_id <- factor(rep(draws, each = 100 * nrow(newdat)))
    df_full$ellipse_id <- factor(rep(1:(n * nrow(newdat)), each = 100))

    dfs_unlist1 <- lapply(
      dfs_list,
      Reduce,
      f = rbind
    )

    return(
      cbind(
        df_full,
        Reduce(rbind, dfs_unlist1)
      )
    )

  }
}


# helper functions ---------------------------

rhs_form <- function(form){
  sform <- Reduce(paste, deparse(form))
  formula(stringr::str_extract(sform, "\\~ .+"))
}


make_L2d <- function(rho, sigma){
  Omeg <- diag(nrow = 2)
  Omeg[upper.tri(Omeg, diag = F)] <- rho
  Omeg[lower.tri(Omeg, diag = F)] <- rho
  S <- diag(sigma) %*% Omeg %*% diag(sigma)
  return(t(chol(S)))
}

make_parmat <- function(P, theta_1, theta_2){

  B <- matrix(data = 0, nrow = sum(P), ncol = 2)
  B[1:P[1], 1] <- theta_1
  B[(P[1] + 1):sum(P), 2] <- theta_2
  return(B)

}

ellipse_df_list <- function(mu_new, L_new, ones, xy_circ){
  mapply(
    function(mu, L, ones, xy_circ){
      t(diag(mu) %*% ones + L %*% xy_circ)
    },
    mu_new, L_new,
    MoreArgs = list(
      ones = ones,
      xy_circ = xy_circ
    ),
    SIMPLIFY = F
  )
}





