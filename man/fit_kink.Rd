\name{fit_kink}
\alias{fit_kink}
\title{Fit a Right-Side Functional Kink Model}
\usage{
fit_kink(data, gamma_grid = NULL, grid_len = 300, nbasis = 10,
  rho = 1e-3, degree = 3, diff_order = 2)
}
\arguments{
\item{data}{A \code{kink_data} object returned by \code{prep_kink_data}.}
\item{gamma_grid}{Optional numeric grid of candidate kink locations.}
\item{grid_len}{Number of grid points when \code{gamma_grid} is not supplied.}
\item{nbasis}{Number of B-spline basis functions.}
\item{rho}{Roughness penalty strength.}
\item{degree}{B-spline degree.}
\item{diff_order}{Difference order used in the roughness penalty.}
}
\description{
Fits the right-side kink model by grid search over \eqn{\gamma} and penalized
least squares for functional coefficient curves.
}
\value{A \code{kink_fit} object containing \code{gamma}, \code{beta}, residuals, fitted values, and the RSS path.}
