\name{test_kink}
\alias{test_kink}
\title{Functional Kink Effect Test}
\usage{
test_kink(data, B = 500, nbasis = 12, grid_len = 100, seed = 123)
}
\arguments{
\item{data}{A \code{kink_data} object returned by \code{prep_kink_data}.}
\item{B}{Number of multiplier bootstrap replications.}
\item{nbasis}{Number of B-spline basis functions under the null model.}
\item{grid_len}{Number of candidate kink locations used in the supremum statistic.}
\item{seed}{Random seed for multiplier bootstrap weights.}
}
\description{
Tests the kink effect using a score/supremum statistic and multiplier bootstrap.
}
\value{A \code{kink_test} object containing \code{p_value}, observed statistic, bootstrap statistics, and the gamma grid.}
