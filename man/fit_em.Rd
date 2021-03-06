% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_algorithm.R
\name{fit_em}
\alias{fit_em}
\title{Fit coding models via an EM algorithm}
\usage{
fit_em(
  long_labels,
  model = "BACE",
  max_iter = 1000,
  tol = 1e-06,
  verbose = FALSE
)
}
\arguments{
\item{long_labels}{should be a data frame with columns:
\itemize{
\item{"ii"}{object index}
\item{"jj"}{coder index}
\item{"yy"}{the coder's label for that object}
}
because long_labels is in long format, there shouldn't be any missing data}

\item{model}{one of c("BACE", "DS", "MACE")}

\item{max_iter}{maximum number of iterations for the EM algorithm; I would advise not touching this unless you know what you are doing.}

\item{tol}{the algorithm is terminated once the change in average log likelihood (average taken at the object level) is both positive and less than tol.}

\item{verbose}{if TRUE, the function will print status updates during the EM algorithm.}
}
\description{
Estimates parameters of either a BACE, MACE, or DS coding model using coder-object labels
}
\keyword{coding}
