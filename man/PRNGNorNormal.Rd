% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRNGNorNormal.R
\name{PRNGNorNormal}
\alias{PRNGNorNormal}
\title{Pseudo-Random Number Generation of Correlated Continuous (Normal-Normal) Outcomes}
\usage{
PRNGNorNormal(M, mvec, mu_r, s_r, s, beta0, beta1, x, seed)
}
\arguments{
\item{M}{The number of groups in the resulting data.}

\item{mvec}{A vector indicating the length of each group in the data.}

\item{mu_r}{The mean of the normal random effect.}

\item{s_r}{The standard deviation of the normal random effect.}

\item{s}{The standard deviation of the random error associated with the response.}

\item{beta0}{The true intercept in the systematic component of the model.}

\item{beta1}{The true slope in the systematic component of the model.}

\item{x}{A vector of predictor values.}

\item{seed}{The seed for data generation.}
}
\description{
This function produces pseudo-random numbers for responses for correlated (longitudinal) data.  The responses are generated as continuous (Normal), using a normal random effect to induce autocorrelation.  The predictors must be provided as a vector.  The function returns a response vector, labeled "Outcomes".
}
\examples{
PRNGNorNormal()
}
\keyword{Correlated}
\keyword{Data}
\keyword{Generation}
\keyword{Simulation}
