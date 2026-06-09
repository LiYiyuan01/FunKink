\name{raw_gamma}
\alias{raw_gamma}
\title{Convert Gamma Back to the Raw X Scale}
\usage{
raw_gamma(gamma, data)
}
\arguments{
\item{gamma}{Estimated gamma on the model scale.}
\item{data}{A \code{kink_data} object.}
}
\description{
Converts an estimated kink location back to the raw x scale when \code{scale_xz = TRUE} was used.
}
\value{Numeric gamma on the raw x scale.}
