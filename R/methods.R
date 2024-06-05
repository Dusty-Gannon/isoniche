

#' Summarize fitted isoniche model
#'
#' @param mfit A fitted isoniche object
#'
#' @return This general method returns a summary of the results of the fitted model, including
#' 95% credible intervals for model coefficients as well as Gelman's Rhat statistic for
#' estimating chain convergence.
#' @export
#'
summary.isoniche <- function(mfit){

  # get group names
  resp_names <- unlist(lapply(
    mfit$model$mean,
    function(x){
      all.vars(x)[1]
    }
  ))

  # first get parameter names as stats package calls them
  mparnames <- lapply(
    mfit$model$mean,
    function(x, dat){
      colnames(
        stats::model.matrix(rhs_form(x), data = dat)
      )
    },
    dat = mfit$data$df
  )

  vparnames <- lapply(
    mfit$model$var,
    function(x, dat){
      colnames(
        stats::model.matrix(x, data = dat)
      )
    },
    dat = mfit$data$df
  )

  joint_post <- rstan::extract(mfit$fit)
  rhat_all <- rstan::summary(mfit$fit)$summary[, "Rhat"]

  # remove lp from list
  joint_post <- joint_post[-length(joint_post)]
  rhat <- rhat_all[-length(rhat_all)]

  # add stats-like names
  for(i in 1:length(mparnames)){
    colnames(joint_post[[i]]) <- mparnames[[i]]
  }
  for(i in 1:length(vparnames)){
    colnames(joint_post[[length(mparnames) + i]]) <- vparnames[[i]]
  }

  # compile into dataframe
  parnames <- c(
    paste("mean", resp_names, sep = "_"),
    paste("sd", resp_names, sep = "_"),
    "corr"
  )
  df_sum <- data.frame(
    coefficient = unlist(lapply(joint_post, colnames)),
    parameter = c(
      rep(parnames, sapply(joint_post, ncol))
    ),
    posterior_mean = unlist(lapply(joint_post, colMeans)),
    q.025 = unlist(lapply(
      joint_post,
      function(x){
        apply(x, 2, quantile, probs = 0.025)
      }
    )),
    q.975 = unlist(lapply(
      joint_post,
      function(x){
        apply(x, 2, quantile, probs = 0.975)
      }
    )),
    Rhat = rhat
  )

  # compile the return list
  isosum <- list(
    model = list(
      model_mean = mfit$model$mean,
      model_var = mfit$model$var[1:2],
      model_corr = mfit$model$var[3]
    ),
    note = "Note: sd_ coefficients are on the log scale and correlation coefficients are on a modified logit scale.\n",
    summary = df_sum
  )

  class(isosum) <- "summary.isoniche"
  isosum
}



#' Predict method for isoniche objects
#'
#' This function can be used to create mean vectors and covariance matrices
#' for groups or conditions specified in \code{newdat}, or create draws from
#' the posterior predictive distribution for conditions specified in \code{newdat}.
#'
#' @param mfit Fitted isoniche model.
#' @param newdat A dataframe specifyng the conditions for which to create the
#' predictions or posterior draws for the mean vector and covariance matrix.
#' @param n Number of draws from the posterior to take. If \code{n = 1}, then
#' the marginal means for each parameter will be used to construct the mean vector
#' and covariance matrix.
#' @param summarize Logical indicating whether to summarize the draws from the posterior or
#' posterior predictive distribution. If \code{TRUE}, the function will return a data frame
#' summarizing the draws.
#' @param type Character of either \code{"mean"} or \code{"response"}. If \code{"mean"}, the
#' function will create draws from the posteriors of the means and variances of each row in
#' \code{newdat}. If \code{"response"}, the function will create \code{npp} draws from the
#' posterior predictive distribution for each of \code{n} draws from the posterior for each
#' row of \code{newdat}.
#' @param npp Number of draws from the posterior predictive distribution for each draw from
#' the joint posterior. These can be used for posterior predictive checks of model assumptions.
#'
#' @importFrom stats model.matrix quantile formula
#' @return Either a dataframe (when \code{summarize = TRUE}) or a list of \code{n} draws from the
#' joint posterior for each row of \code{newdat}.
#' @export
#'
predict.isoniche <- function(mfit, newdat, n = 250, summarize = FALSE, type = "mean", npp = 1){
  # first check names
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

  # construct parameters
  if(n == 1){
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

    # make means and covariances for each row of newdat
    mu_new <- lapply(
      1:nrow(newdat),
      function(i, X_new, B, resp_names){
        mu <- as.double(X_new[i, ] %*% B)
        names(mu) <- resp_names
        return(mu)
      },
      X_new, B, resp_names
    )
    Sigma_new <- lapply(
      1:nrow(newdat),
      function(i, Z_new, Zeta, G_new, gamma, resp_names){
        s <- exp(as.double(Z_new[i, ] %*% Zeta))
        rho <- 2 * stats::plogis(as.double(G_new[i, ] %*% gamma)) - 1
        Sigma <- diag(s^2)
        Sigma[1, 2] <- rho * prod(s)
        Sigma[2, 1] <- rho * prod(s)
        rownames(Sigma) <- resp_names
        colnames(Sigma) <- resp_names
        return(Sigma)
      },
      Z_new, Zeta, G_new, gamma, resp_names
    )
  }

  # now the case with multiple draws from the posterior
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
    # covariance matrix
    mu_new <- lapply(
      B_list,
      function(B, X_new){
        lapply(
          1:nrow(X_new),
          function(i, X_new, B, resp_names){
            mu <- as.double(X_new[i, ] %*% B)
            names(mu) <- resp_names
            return(mu)
          },
          X_new = X_new, B = B, resp_names
        )
      },
      X_new = X_new
    )

    Sigma_new <- mapply(
      function(Z, Zeta, G, gamma, newdat){
        lapply(
          1:nrow(newdat),
          function(i, Z, Zeta, G, gamma, resp_names){
            s <- exp(as.double(Z[i, ] %*% Zeta))
            rho <- 2 * stats::plogis(as.double(G[i, ] %*% gamma)) - 1
            Sigma <- diag(s^2)
            Sigma[1, 2] <- rho * prod(s)
            Sigma[2, 1] <- rho * prod(s)
            rownames(Sigma) <- resp_names
            colnames(Sigma) <- resp_names
            return(Sigma)
          },
          Z, Zeta, G, gamma, resp_names
        )
      },
      Zeta_list, gamma_list,
      MoreArgs = list(
        Z = Z_new, G = G_new, newdat = newdat
      ),
      SIMPLIFY = F
    )
  }

  # return predictions
  if(type == "mean" & n == 1){
    res <- lapply(
      1:nrow(newdat),
      function(i, mu, Sigma, newdat){
        return(list(
          cond = newdat[i, ],
          mu = as.double(mu[[i]]),
          Sigma = Sigma[[i]]
        ))
      },
      mu_new, Sigma_new, newdat
    )
  }

  if(type == "mean" & n > 1){
    # compile into output
    res <- lapply(
      1:nrow(newdat),
      function(i, mu, Sigma, newdat){
        retlist <- list(
          cond = newdat[i, ],
          mu = t(sapply(mu, function(x, i){x[[i]]}, i = i)),
          Sigma = lapply(Sigma, function(x, i){x[[i]]}, i = i)
        )
      },
      mu_new, Sigma_new, newdat
    )
    if(summarize){
      res_sum1 <- lapply(
        res,
        function(x){
          mu_sum <- data.frame(
            param = paste("mu", colnames(x$mu), sep = "_"),
            mean = colMeans(x$mu),
            q025 = apply(x$mu, 2, stats::quantile, probs = 0.025),
            q975 = apply(x$mu, 2, stats::quantile, probs = 0.975)
          )
          S_flat <- t(sapply(x$Sigma, as.vector))[, -2]
          S_sum <- data.frame(
            param = c(
              paste0("var_", colnames(x$mu)[1]),
              paste(c("cov", colnames(x$mu)), collapse = "_"),
              paste0("var_", colnames(x$mu)[2])
            ),
            mean = colMeans(S_flat),
            q025 = apply(S_flat, 2, stats::quantile, probs = 0.025),
            q975 = apply(S_flat, 2, stats::quantile, probs = 0.975)
          )
          return(rbind(mu_sum, S_sum))
        }
      )
      res_sum2 <- lapply(
        1:length(res),
        function(i, res_sum1){
          res[[i]]$cond[rep(1, nrow(res_sum1[[i]])), ]
        },
        res_sum1
      )
      res <- cbind(
        Reduce(rbind, res_sum2),
        Reduce(rbind, res_sum1)
      )
    }
  }

  # now do the case for posterior predictions
  if(type == "response"){
    if(n < 10){
      warning("n may be too small to capture the variability in the posterior predictive
              distribution. Consider increasing.")
    }
    y_new <- lapply(
      1:nrow(newdat),
      function(i, mu_new, Sigma_new, n){
        y_new <- mapply(
          function(mu_new, Sigma_new, n, i){
            mvtnorm::rmvnorm(n, mean = mu_new[[i]], sigma = Sigma_new[[i]])
          },
          mu_new, Sigma_new, n, i,
          SIMPLIFY = F
        )
      },
      mu_new, Sigma_new, n
    )
    y_new <- lapply(
      y_new,
      function(x, names){
        mat <- Reduce(rbind, x)
        colnames(mat) <- names
        return(mat)
      },
      names = resp_names
    )

    if(summarize){
      y_new <- lapply(
        y_new,
        function(x){
          df <- data.frame(
            var = colnames(x),
            mean = colMeans(x),
            q.025 = apply(x, 2, stats::quantile, probs = 0.025),
            q.1 = apply(x, 2, stats::quantile, probs = 0.1),
            q.9 = apply(x, 2, stats::quantile, probs = 0.9),
            q.975 = apply(x, 2, stats::quantile, probs = 0.975)
          )
          return(df)
        }
      )
    }

    # compile results
    res <- lapply(
      1:nrow(newdat),
      function(i, y_new, newdat){
        df <- newdat[rep(i, nrow(y_new[[i]])), ]
        return(cbind(df, y_new[[i]]))
      },
      y_new, newdat
    )
  }

  return(res)

}




#' Print methods for the summarized isoniche object
#'
#' @param sumobj summary.isoniche object
#' @param digits Number of digits to include in printed table. Defaults to 3.
#'
#' @export
#'
print.summary.isoniche <- function(sumobj, digits = 3){
    df_sum_print <- sapply(
      1:ncol(sumobj$summary),
      function(i, df){
        if(is.double(df[,i])){
          return(round(df[,i], digits))
        } else{
          return(df[,i])
        }
      },
      df = sumobj$summary
    )
    df_sum_print <- as.data.frame(df_sum_print)
    names(df_sum_print) <- names(sumobj$summary)

    # first mean model
    cat(
      c(
        "Model for mean vector:",
        as.character(unlist(sumobj$model$model_mean)),
        "link: identity\n"
      ),
      sep = "\n"
    )

    # print sd model
    cat(
      c(
        "Model for scale (SD) vector:",
        as.character(unlist(sumobj$model$model_var)),
        "link: log()\n"
      ),
      sep = "\n"
    )

    # print model for correlation
    cat(
      c(
        "Model for correlation:",
        as.character(unlist(sumobj$model$model_corr)),
        "link: logit(0.5 * (sigma + 1))\n"
      ),
      sep = "\n"
    )

    print(df_sum_print)

}


#' Mean of posterior normalized residual distribution
#'
#' @param mfit Fitted isoniche object
#'
#' @return Matrix of posterior mean residuals, one column for each response
#' variable.
#' @export
#'
residuals.isoniche <- function(mfit){
  # posteriors of fitted values
  y_hat <- predict(mfit, newdat = mfit$data$df, n = 100)

  mu <- lapply(
    y_hat,
    function(x){return(x$mu)}
  )

  # put y vecs into a list
  y <- lapply(
    1:nrow(mfit$data$datlist$y),
    function(i, ymat){
      return(as.double(ymat[i, ]))
    },
    ymat = mfit$data$datlist$y
  )

  # extract names of response variables
  resp_names <- unlist(lapply(
    mfit$model$mean,
    function(x){
      all.vars(x)[1]
    }
  ))

  # compute inverse square root of covariance matrix
  Ls <- lapply(
    y_hat,
    function(x){
      lapply(
        x$Sigma,
        function(S){ solve(t(chol(S))) }
      )
    }
  )

  # get and normalize residuals
  resids <- t(mapply(
    function(y, mu, L_list){
      res_i <- lapply(
        1:length(L_list),
        function(i, y, mu, L_list){
          as.double(L_list[[i]] %*% (y - mu[i, ]))
        },
        y = y, mu = mu, L_list = L_list
      )
      # now summarize
      res_i_hat <- colMeans(
        Reduce(rbind, res_i)
      )
      return(res_i_hat)
    },
    y, mu, Ls
  ))

  colnames(resids) <- paste("res", resp_names, sep = "_")

  return(resids)
}


#' Plot posterior means of residual distributions
#'
#' @param mfit Fitted isoniche object
#'
#' @return Multipanel plot of residuals. The first shows the bivariate residuals, which should
#' be approximately spherical (i.e., be bivariate normal with mean zero and covariance matrix
#' \eqn{I_2}). The next plots show qq plots for each dimension.
#' @export
#'
plot.isoniche <- function(mfit){

  resids <- residuals(mfit)

  sd1 <- cbind(
    cos(seq(0, 2 * pi, length.out = 100)),
    sin(seq(0, 2 * pi, length.out = 100))
  )

  sd2 <- 2 * sd1

  par(mfrow = c(2,2), mar = c(4,4,1,1))
  plot(resids, xlab = colnames(resids)[1], ylab = colnames(resids)[2])
  lines(sd1, col = "blue", lty = "dashed")
  lines(sd2, col = "red", lty = "dashed")
  qqnorm(resids[, 1], main = colnames(resids)[1])
  abline(0, 1, col = "blue")
  qqnorm(resids[, 2], main = colnames(resids)[2])
  abline(0, 1, col = "blue")

}


