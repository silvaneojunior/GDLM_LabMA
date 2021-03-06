% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{poisson_pred}
\alias{poisson_pred}
\title{poisson_pred}
\usage{
poisson_pred(model, IC_prob = 0.95)
}
\arguments{
\item{model}{List: contaning the parameters of the conjugated distribuition for the linear predictor (may be based on the prior, filtered or smoothed distribuition).}

\item{IC_prob}{Numeric: the desired credibility for the credibility interval}
}
\value{
A list containing the following values:
\itemize{
   \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
   \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
   \item icl.pred vector/matrix: the percentile of 100*((1-IC_prob)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
   \item icu.pred vector/matrix: the percentile of 100*(1-(1-IC_prob)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
}
}
\description{
Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
It's worth noting that, since the conjugated distribuition of the linear predictior is Gamma (see Ref. Raíra), then the predictive distribuition is Negative Binomial.
}
\examples{
# A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.

model=list(
'a'=c(1:3),
'b'=c(3:1))
)

poisson_pred(model)
}
