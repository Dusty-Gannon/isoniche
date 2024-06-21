
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
  ynames <- sapply(mean, function(x){
    all.vars(x)[1]
  })

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

  # compile standard list for return
  ret_obj <- list(
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
  class(ret_obj) <- "isoniche"
  return(
    ret_obj
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
#' @param sds Multiple of standard deviations to use for ellipses. For standard ellipses, \code{sds = 1},
#' For approximate 95% ellipses, you could use \code{sds = 2}.
#'
#' @return A dataframe that can be used for plotting standard ellipses or calculating isotopic niche
#' statistics.
#'
#' @importFrom stats formula model.matrix
#' @importFrom mvtnorm rmvnorm
#' @export
#' @example /inst/examples/isoniche_eg.R
#'
construct_ellipses <- function(mfit, newdat, n = 1, q = 1){

  resp_names <- unlist(lapply(
    mfit$model$mean,
    function(x){
      all.vars(x)[1]
    }
  ))
  varnames <- unlist(lapply(
    mfit$model$mean,
    function(x){
      all.vars(x)[-1]
    }
  ))
  varnames <- c(
    varnames,
    unlist(lapply(mfit$model$var, all.vars))
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
  ) * sqrt(q)
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


#' Construct posterior of SEA
#'
#' This function takes in a dataframe that represents the groups or conditions over which to compute
#' the standard ellipse area (SEA) and takes \code{n} draws from the posterior covariance matrix
#' for each row of \code{newdat}. It then computes the SEA for each and returns a dataframe
#' of the posterior draws with one row per draw per condition.
#'
#' @param mfit Fitted isoniche model object.
#' @param newdat New data, with all the same columns as the original data used to
#' fit the model.
#' @param n Number of draws from the posterior SEA for each group or condition in \code{newdat}.
#'
#' @return A dataframe that can be used for plotting posterior SEAs or computing summary statistics.
#'
#' @export
#' @example /inst/examples/isoniche_eg.R
sea <- function(mfit, newdat, n = 250){

  preds <- predict.isoniche(mfit, newdat, n)

  draws_sea <- lapply(
    preds,
    function(x){
      sapply(
        x$Sigma,
        FUN = function(s){
          lambda <- eigen(s, only.values = TRUE)
          pi * prod(sqrt(lambda$values))
        }
      )
    }
  )

  # convert to dataframe
  return_df <- newdat[rep(1:nrow(newdat), each = n), ]
  return_df$sea <- unlist(draws_sea)
  return_df$draw_id <- paste("s", rep(1:n, nrow(newdat)), sep = "_")
  return(return_df)

}





#' Calculate the overlap of ellipses using Monte Carlo simulation
#'
#' This function estimates the overlap area of multiple ellipses
#' defined by their means and covariance matrices using a Monte Carlo
#' simulation approach.
#'
#' @param mus A list of numeric vectors representing the mean vectors of the ellipses.
#' @param Sigmas A list of matrices representing the covariance matrices of the ellipses.
#' @param q A numeric value representing the quantile used to define the ellipses.
#' This can be defined based on the desired quantile of a \eqn{\chi^2_2} distribution.
#' @param npp An integer specifying the number of points to be used in the Monte Carlo simulation.
#' Default is 1000.
#'
#' @return A numeric value representing the estimated overlap area of the ellipses.
#'
#' @details
#' The function first simulates `npp` points uniformly within a bounding box that contains
#' the ellipses. It then tests whether these points lie within each of the ellipses defined
#' by their respective means and covariance matrices. The overlap area is calculated based
#' on the proportion of points that lie within all ellipses, multiplied by the area of the
#' bounding box.
#'
#' @examples
#' \dontrun{
#' mus <- list(c(0, 0), c(1, 1))
#' Sigmas <- list(diag(2), diag(2))
#' q <- 2
#' ellipse_overlap_single_mc(mus, Sigmas, q)
#' }
#'
#' @importFrom stats runif
#' @export
ellipse_overlap_single_mc <- function(mus, Sigmas, q, npp = 1000){

  # first, simulate some points in the sample space
  bb <- find_bound_box(mus, Sigmas, q)
  y_sim <- cbind(
    stats::runif(npp, min = bb[1,1], max = bb[2,1]),
    stats::runif(npp, min = bb[2,2], max = bb[3,2])
  )
  ones <- matrix(1, ncol = ncol(y_sim), nrow = nrow(y_sim))

  # get tests
  tests <- mapply(
    FUN = function(mu, Sigma, ones, q, y_sim){
      test <- diag((y_sim - ones %*% diag(mu)) %*% solve(Sigma) %*% (t(y_sim) - diag(mu) %*% t(ones)))
      return(test < q)
    },
    mus, Sigmas, MoreArgs = list(ones = ones, q = q, y_sim = y_sim),
    SIMPLIFY = F
  )

  # now find which sims are in all
  sims_in_all <- Reduce("&", tests)

  area_S <- (bb[2, 1] - bb[1, 1]) * (bb[3, 2] - bb[2, 2])
  return(area_S * mean(sims_in_all))
}





#' Estimate Posterior Overlap of Pairs of Ellipses
#'
#' This function estimates the overlap area for each pair of ellipses
#' defined by the user for each of `n` draws from the posterior distribution
#' of the mean vector and covariance matrix defining each ellipse.
#'
#' @param mfit An `isoniche` model object.
#' @param df_pairs A data frame or a list of data frames. If a list of dataframes,
#'  each should have 2 rows defining the comparisons. Each row represents a condition
#'  for which the ellipse is defined. If a single data frame, the function will compute all pairwise
#'  overlaps by default.
#' @param q The quantile used to define the ellipses. This can be defined based on
#' the desired quantile of a \eqn{\chi^2_2} distribution. The default is 1.
#' @param n An integer specifying the number of draws from the posterior distribution. Default is 250.
#' @param npp An integer specifying the number of points to be used in the Monte Carlo simulation for
#' estimating ellipse overlap. Default is 1000. Overlap estimates will be more precise with more draws,
#' but the function will get slow.
#'
#' @return A list where each element corresponds to a pair of conditions and contains the following components:
#'   - `cond1`: The first condition in the pair.
#'   - `cond2`: The second condition in the pair.
#'   - `overlap_area`: A numeric vector of estimated overlap areas for each draw from the posterior.
#'   This approximates the posterior distribution of the area of overlap.
#'
#' @details
#' For pair of conditions and each draw from the posterior distribution, the function first
#' simulates `npp` points uniformly within a bounding box that contains the ellipses.
#' It then tests whether these points lie within each of the ellipses defined
#' by their respective means and covariance matrices. The overlap area is calculated based
#' on the proportion of points that lie within all ellipses, multiplied by the area of the
#' bounding box.
#'
#'
#' @seealso \code{\link{ellipse_overlap_single_mc}}
#' @export
ellipse_overlap_mc <- function(mfit, df_pairs, q = 1, n = 250, npp = 1000){

  if(!(is.data.frame(df_pairs) | is.list(df_pairs))){
    stop("Unrecognized input type for df_pairs.
         Please provide a data.frame or a list of data.frames each with 2 rows defining the comparisons.\n")
  }
  # convert dataframe into list of pairs
  if(is.data.frame(df_pairs)){
    n_conds <- nrow(df_pairs)
    pairs_list <- vector(mode = "list")
    for(i in 1:(n_conds - 1)){
      for(j in (i + 1):n_conds){
        pairs_list <- c(
          pairs_list,
          list(df_pairs[c(i, j), ])
        )
      }
    }
  } else{
    pairs_list <- df_pairs
  }

  # get draws from posterior means and covariances
  preds_pairs <- lapply(
    pairs_list,
    function(df, mfit, n){
      predict(mfit, df, n, type = "mean")
    },
    mfit, n
  )

  # compute overlap for each of the draws from the posterior
  # and for each pair
  overlap <- lapply(
    preds_pairs,
    FUN = function(pair, q, npp){
      lapply(
        1:n,
        FUN = function(i, l1, l2, q, npp){
          ellipse_overlap_single_mc(
            mus = list(l1$mu[i, ], l2$mu[i, ]),
            Sigmas = list(l1$Sigma[[i]], l2$Sigma[[i]]),
            q = q,
            npp = npp
          )
        },
        l1 = pair[[1]], l2 = pair[[2]],
        q = q, npp = npp
      )
    },
    q = q, npp = npp
  )

  # add dataframe listing the comparison
  for(i in 1:length(pairs_list)){
    overlap[[i]] <- c(
      list(
        cond1 = pairs_list[[i]][1, ], cond2 = pairs_list[[i]][2, ],
        overlap_area = unlist(overlap[[i]])
      )
    )
  }

  return(overlap)
}


# ellipse_overlap <- function(mu1, mu2, Sigma1, Sigma2, q = 1){
#
#   # first invert covariance matrices
#   A1 <- solve(Sigma1)
#   A2 <- solve(Sigma2)
#
#   b1 <- -mu1
#   b2 <- -mu2
#
#   # next find intersections of ellipses
#
#   # define quartic function in y
#   quart_in_y <- function(x, u){
#     u[1] + u[2] * x + u[3] * x^2 + u[4] * x^3 + u[5] * x^4
#   }
#
#   # defining v
#   v <- vector("double", length = 11)
#   v[1] <- 2 * (A1[1,1] * A2[1,2] - A2[1,1] * A1[1,2])
#   v[2] <- A1[1,1] * A2[2,2] - A2[1,1] * A1[2,2]
#   v[3] <- A1[1,1] * b2[1] - A2[1,1] * b1[1]
#   v[4] <- A1[1,1] * b2[2] - A2[1,1] * b1[2]
#   v[5] <- q * (A1[1,1] - A2[1,1])
#   v[6] <- 2 * (A1[1,2] * A2[2,2] - A2[1,2] * A1[2,2])
#   v[7] <- 2 * (A1[1,2] * b2[2] - A2[1,2] * b1[2])
#   v[8] <- 2 * q * (A1[1,2] - A2[1, 2])
#   v[9] <- A1[2,2] * b2[1] - A2[2,2] * b1[2]
#   v[10] <- b1[1] * b2[2] - b2[1] * b1[2]
#   v[11] <- q * (b1[1] - b2[1])
#
#   # now define u based on v
#   u <- double(5)
#   u[1] <- v[3] * v[11] - v[5]^2
#   u[2] <- v[1] * v[11] + v[3] * (v[8] + v[10]) - 2 * v[4] * v[5]
#   u[3] <- v[1] * (v[8] + v[10]) + v[3] * (v[7] - v[9]) - v[4]^2 - 2 * v[2] * v[5]
#   u[4] <- v[1] * (v[7] - v[9]) + v[3] * v[6] - 2 * v[2] * v[4]
#   u[5] <- v[1] * v[6] - v[2]^2
#
#   # now find roots of the quartic
#   y_star <- uniroot(quart_in_y, interval = c(-2, 2), u = u)
#
#
# }


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

find_bound_box <- function(mus, Sigmas, q){
  xy_circ <- rbind(
    cos(seq(0, 2 * pi, length.out = 100)),
    sin(seq(0, 2 * pi, length.out = 100))
  ) * sqrt(q)
  Ls <- lapply(Sigmas, function(x){t(chol(x))})
  ones <- matrix(data = 1, 2, 100)

  dfs_list <- ellipse_df_list(mus, Ls, ones, xy_circ)

  # now get the min and max of each coordinate
  all <- Reduce(rbind, dfs_list)

  bb <- rbind(
    c(min(all[, 1]), min(all[, 2])),
    c(max(all[, 1]), min(all[, 2])),
    c(max(all[, 1]), max(all[, 2])),
    c(min(all[, 1]), max(all[, 2]))
  )

  bb[c(1,4), 1] <- bb[c(1,4), 1] - (1/20) * (bb[2,1] - bb[1,1])
  bb[c(2,3), 1] <- bb[c(2,3), 1] + (1/20) * (bb[2,1] - bb[1,1])
  bb[c(1,2), 2] <- bb[c(1,2), 2] - (1/20) * (bb[4,2] - bb[1,2])
  bb[c(3,4), 2] <- bb[c(3,4), 2] + (1/20) * (bb[4,2] - bb[1,2])

  return(bb)
}

