

summary.isoniche <- function(mfit){

  # get the call
  call <- deparse(substitute(summary(mfit)))
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

  # remove lp from list
  joint_post <- joint_post[-length(joint_post)]

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
    "cor"
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
    ))
  )

  if(interactive()){
    df_sum_print <- sapply(
      1:ncol(df_sum),
      function(i, df){
        if(is.double(df[,i])){
          return(round(df[,i], 2))
        } else{
          return(df[,i])
        }
      },
      df = df_sum
    )
    df_sum_print <- as.data.frame(df_sum_print)
    names(df_sum_print) <- names(df_sum)
    print(df_sum_print)
    cat("Note: sd_ coefficients are on the log scale and correlation coefficients are on a modified logit scale.\n")
  }

  invisible(df_sum)

}
