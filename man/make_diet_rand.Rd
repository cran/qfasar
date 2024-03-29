% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_diet_rand.R
\name{make_diet_rand}
\alias{make_diet_rand}
\title{Generate random diet compositions}
\usage{
make_diet_rand(uniq_types, n_diet)
}
\arguments{
\item{uniq_types}{A factor of unique prey-type names.}

\item{n_diet}{The integer number of diet compositions to generate.}
}
\value{
A list containing the following elements: \describe{
  \item{diet_rand}{A numeric matrix of random diet compositions, in
    column-major format.}
  \item{err_code}{An integer error code (0 if no error is detected).}
  \item{err_message}{A string contains a brief summary of the execution.}
  }
}
\description{
The function \code{make_diet_rand} generates a user-specified number of
random diet compositions.
}
\section{Details}{

The function \code{make_diet_rand} generates a specified number of random
diet compositions to support simulation-based research of the performance
of QFASA diet estimation procedures. Given a diet composition, predator
fatty acid signatures can be generated using the
function \code{\link{make_pred_sigs}}. The diets of such simulated predators
can then be estimated, and the diet estimates can be compared to the
known diet composition to evaluate bias, variance, and perhaps other
properties.

The algorithm starts by generating a uniformly distributed random number
between 0 and 1 as the diet proportion for the first prey type. The algorithm
then considers each additional prey type in turn, generating a uniform random
number between zero and 1 minus the sum of the proportions assigned to the
preceding prey types. The diet proportion for the last prey type is 1 minus
the sum of the other diet proportions. As a hedge against limitations in the
random number generator, the proportions are then randomly ordered among prey
types.

It is critical that the prey-type names match those in the prey library. The
easiest way to ensure this happens is to pass the object uniq_types returned
a call to the function \code{\link{prep_sig}} as the uniq_types argument.
Alternatively, and more risky, a vector of unique prey names can be created
using the concatenate function and cast as a factor, i.e.,
uniq_types <- as.factor(c("Prey_1", "Prey_2", ..., "Prey_P)).
}

\examples{
make_diet_rand(uniq_types = as.factor(c("Bearded",
                                        "Beluga",
                                        "Bowhead",
                                        "Ribbon",
                                        "Ringed",
                                        "Spotted",
                                        "Walrus")),
               n_diet = 100)

}
