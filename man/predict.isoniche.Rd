% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{predict.isoniche}
\alias{predict.isoniche}
\title{Predict method for isoniche objects}
\usage{
\method{predict}{isoniche}(mfit, newdat, n = 250, summarize = FALSE, type = "mean", npp = 1)
}
\arguments{
\item{mfit}{Fitted isoniche model.}

\item{newdat}{A dataframe specifyng the conditions for which to create the
predictions or posterior draws for the mean vector and covariance matrix.}

\item{n}{Number of draws from the posterior to take. If \code{n = 1}, then
the marginal means for each parameter will be used to construct the mean vector
and covariance matrix.}

\item{summarize}{Logical indicating whether to summarize the draws from the posterior or
posterior predictive distribution. If \code{TRUE}, the function will return a data frame
summarizing the draws.}

\item{type}{Character of either \code{"mean"} or \code{"response"}. If \code{"mean"}, the
function will create draws from the posteriors of the means and variances of each row in
\code{newdat}. If \code{"response"}, the function will create \code{npp} draws from the
posterior predictive distribution for each of \code{n} draws from the posterior for each
row of \code{newdat}.}

\item{npp}{Number of draws from the posterior predictive distribution for each draw from
the joint posterior. These can be used for posterior predictive checks of model assumptions.}
}
\value{
Either a dataframe (when \code{summarize = TRUE}) or a list of \code{n} draws from the
joint posterior for each row of \code{newdat}.
}
\description{
This function can be used to create mean vectors and covariance matrices
for groups or conditions specified in \code{newdat}, or create draws from
the posterior predictive distribution for conditions specified in \code{newdat}.
}
