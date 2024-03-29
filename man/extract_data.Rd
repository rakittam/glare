% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{extract_data}
\alias{extract_data}
\title{Data extraction}
\usage{
extract_data(object, newdata = NULL, parameter = NULL, ...)
}
\arguments{
\item{object}{an object of class \code{"\link{glare}"}.}

\item{newdata}{an optional data frame, list or environment containing the
variables in the model. If not used, the variables are taken from the
objective.}

\item{parameter}{optional parameter input. Default takes the parameter
from the glare objective.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
list of data extracted from input.
}
\description{
Returns data from object input.
}
