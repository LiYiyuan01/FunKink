\name{prep_kink_data}
\alias{prep_kink_data}
\title{Prepare Data for a Right-Side Functional Kink Model}
\usage{
prep_kink_data(Y, x, Z = NULL, t = NULL, scale_xz = TRUE,
  center_y = FALSE, standardize_y = FALSE, intercept = FALSE,
  remove_na = TRUE)
}
\arguments{
\item{Y}{Numeric matrix of functional responses, with rows indexing subjects and columns indexing grid points.}
\item{x}{Numeric vector of scalar covariates.}
\item{Z}{Optional numeric matrix of additional scalar covariates. Include a column of ones here if following the simulation convention.}
\item{t}{Optional grid for the functional domain. Defaults to an equally spaced grid on [0,1].}
\item{scale_xz}{Whether to standardize x and Z.}
\item{center_y}{Whether to center Y columnwise.}
\item{standardize_y}{Whether to standardize Y columnwise.}
\item{intercept}{Whether to add an intercept column internally. Use FALSE if Z already contains a column of ones.}
\item{remove_na}{Whether to remove rows with missing values.}
}
\description{
Prepares data for the model
\deqn{Y_i(t) = \beta_1(t)X_i + \beta_2(t)(X_i-\gamma)_+ + Z_i^\top\beta_3(t) + e_i(t).}
}
\value{A list of class \code{kink_data}.}
