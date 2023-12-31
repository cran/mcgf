% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcgf_rs.R
\name{new_mcgf_rs}
\alias{new_mcgf_rs}
\title{Create an mcgf_rs object}
\usage{
new_mcgf_rs(x, label)
}
\arguments{
\item{x}{An mcgf object.}

\item{label}{A vector of regime labels. Its length must be the same as
the number rows in \code{data}.}
}
\value{
An S3 object of class \code{mcgf_rs}. As it inherits and extends the
\code{mcgf} and then the\code{data.frame} class, all methods remain valid to the
\code{data} part of the object. Additional attributes may be assigned and
extracted.
}
\description{
Create an mcgf_rs object
}
\keyword{internal}
