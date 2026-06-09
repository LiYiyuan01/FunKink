\name{predict_kink}
\alias{predict_kink}
\title{Predict from a Functional Kink Fit}
\usage{
predict_kink(object, newdata = NULL, ...)
}
\arguments{
\item{object}{A \code{kink_fit} object.}
\item{newdata}{Optional list with elements \code{x} and \code{Z}. If omitted, fitted values are returned.}
\item{...}{Reserved for future use.}
}
\description{
Returns fitted values for training data or predictions for new data.
}
\value{A numeric matrix of predicted functional responses.}
