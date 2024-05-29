

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


#
# predict.isoniche <- function(mfit, newdat, summarize = FALSE){
#
# }
#



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




