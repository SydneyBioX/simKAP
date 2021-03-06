% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state_algorithm0.R
\name{state_balance}
\alias{state_balance}
\title{state balance function}
\usage{
state_balance(SimResults, recip_matrix = NULL, donor_matrix = NULL)
}
\arguments{
\item{SimResults}{a data.frame, this is the current simulation results}

\item{recip_matrix}{a data.frame, this is the recipient matrix}

\item{donor_matrix}{a data.frame, this is the donor vector}
}
\value{
a vector, which contains all the kidneys allocated to another state
}
\description{
this is the function to obtain the number of kidneys allocated to another state
}
\examples{
data("rawdata", package = "simKAP")
data("newdata", package = "simKAP")

state_balance <- state_balance(rawdata,raw_recip_matrix_subset,raw_donor_matrix[1,]);
}
