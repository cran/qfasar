% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diet_pool.R
\name{diet_pool}
\alias{diet_pool}
\title{Pool diet estimates to combined prey types}
\usage{
diet_pool(rep_grp, est_ind, var_ind = NA, est_mean = NA, var_mean = NA)
}
\arguments{
\item{rep_grp}{The post-multiplication matrix returned by a call to
\code{\link{make_prey_part}} as the object \code{pool_post}, or a
user-defined matrix for custom pooling. Each column defines a prey type
to which estimates should be pooled.}

\item{est_ind}{A numeric matrix of the estimated diet compositions of
individual predators using a partitioned prey library, intended to be the
object est_ind returned by a call to \code{\link{est_diet}}.}

\item{var_ind}{A numeric array containing the estimated variance matrix for
the estimated diet of each predator, intended to be the object
var_ind returned by a call to \code{\link{est_diet}}. Optional.}

\item{est_mean}{A numeric matrix containing the estimated mean diet of each
predator type, intended to be the object est_mean returned by a call to
\code{\link{est_diet}}. Optional.}

\item{var_mean}{A numeric array containing the estimated variance matrix for
the estimated mean diet of each predator type, intended to be the object
var_mean returned by a call to \code{\link{est_diet}}. Optional.}
}
\value{
A list containing the following elements, all of which are organized
  on the basis of the original, unpartitioned prey types: \describe{
  \item{est_ind}{A numeric matrix of the estimated diet compositions of
  individual predators.}
  \item{var_ind}{A numeric array containing the estimated variance matrix for
  the estimated mean diet of each predator.}
  \item{est_mean}{A numeric matrix containing the estimated mean diet of each
    predator type.}
  \item{var_mean}{A numeric array containing the estimated variance matrix
    for the estimated mean diet of each predator type.}
  \item{err_code}{An integer error code (0 if no error is detected).}
  \item{err_message}{A string containing a brief summary of the results.}
}
}
\description{
\code{diet_pool} pools estimated diets and variance matrices to a smaller
number of combined prey types. If \code{\link{est_diet}} is used to estimate
predator diet composition using a partitioned prey library
(\code{\link{make_prey_part}}), \code{diet_pool} pools the partitioned
results back to the original, unpartitioned prey types.
}
\section{Details}{

The function \code{\link{dimac}} explores the prey library for additional
structure with identified prey types. If significant structure is found
within a library, estimating diet composition on the basis of a
partitioned prey library may lead to estimates with less bias and possibly
less variation through reduced prey confounding (Bromaghin et al. 2016). The
function \code{\link{make_prey_part}} takes the clustering results returned
by \code{\link{dimac}} and user specification of the number of clusters in
which to partition each prey type and returns a partitioned prey library
that is ready for use in diet estimation.

However, when estimating diet composition using a partitioned prey library
one may still wish to pool partitioned estimates back to the original,
unpartitioned prey types for reporting purposes. That is the purpose of the
function \code{diet_pool}.

NOTE: \code{diet_pool} can also be used to pool estimates into a smaller
number of combined prey types for reporting purposes. For example, imagine a
prey library with a large number of prey types.  If subsets of the prey
types have similar ecological function, their signatures may share some
similarities (prey confounding, Bromaghin et al. 2016).  In such a case, one
may wish to estimate diet on the basis of the full prey library, but
subsequently pool the resulting estimates to a smaller number of combined
prey types for reporting purposes (reporting groups, Bromaghin 2008) to
reduce the effect of prey confounding. \code{diet_pool} can also be used for
this purpose, though the user would need to manually construct the reporting
group matrix \code{rep_grp}.
}

\section{References}{

Bromaghin, J.F. 2008. BELS: Backward elimination locus selection for studies
  of mixture composition or individual assignment. \emph{Molecular Ecology
  Resources} 8:568-571.

Bromaghin, J.F., S.M. Budge, and G.W. Thiemann. 2016. Should fatty
  acid signature proportions sum to 1 for diet estimation?
  \emph{Ecological Research} 31:597-606.
}

\examples{
diet_pool(rep_grp = matrix(c(1, 0, 0, 0, 0, 0, 0,
                             0, 1, 0, 0, 0, 0, 0,
                             0, 1, 0, 0, 0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0,
                             0, 0, 0, 1, 0, 0, 0,
                             0, 0, 0, 1, 0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0,
                             0, 0, 0, 0, 0, 1, 0,
                             0, 0, 0, 0, 0, 1, 0,
                             0, 0, 0, 0, 0, 0, 1),
                           nrow = 10, byrow = TRUE),
         est_ind = matrix(c(0.116, 0.315,
                            0.028, 0.073,
                            0.000, 0.000,
                            0.131, 0.120,
                            0.000, 0.000,
                            0.000, 0.000,
                            0.723, 0.452,
                            0.000, 0.000,
                            0.000, 0.000,
                            0.002, 0.040),
                          nrow = 10, byrow = TRUE))


}
