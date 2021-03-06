% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_sigs_2_mean.R
\name{dist_sigs_2_mean}
\alias{dist_sigs_2_mean}
\title{Distance between fatty acid signatures and their mean}
\usage{
dist_sigs_2_mean(sig_data, dist_meas = 1, gamma = 1)
}
\arguments{
\item{sig_data}{A numeric matrix of fatty acid signatures in column-major
format.}

\item{dist_meas}{An integer indicator of the distance measure to compute.
Default value 1.}

\item{gamma}{The power parameter of the chi-square distance measure. Default
value 1.}
}
\value{
A list containing the following elements:\describe{
  \item{dist_sum}{The summed distance between each signature and the mean
    signature.}
  \item{dist_mean}{The mean distance between each signature and the mean
    signature.}
  }
}
\description{
The utility function \code{dist_sigs_2_mean} computes the summed and mean
distance between a collection of fatty acid signatures and their mean
signature.
}
\section{Details}{

This is an internal utility function.  The signature data in \code{sig_data}
are presumed to be ready for analysis, which is best accomplished by a call
to the function \code{\link{prep_sig}}.  Consequently, to increase execution
speed during simulations, no numeric error checking is performed.  Please
refer to the documentation for \code{\link{prep_sig}} for additional
information.

Please refer to the documentation for \code{\link{dist_between_2_sigs}} for
additional information regarding distance measures.

Utility functions called by \code{dist_sigs_2_mean}:
  \itemize{
  \item \code{\link{dist_between_2_sigs}}
  }
}

