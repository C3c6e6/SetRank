\name{createIDConverter}
\alias{createIDConverter}
\title{Create a function to convert gene or protein IDs}
\usage{
createIDConverter(annotationPackageName, from, to)
}
\arguments{
  \item{annotationPackageName}{The name of the
  \code{AnnotationDb} package of that will be used to
  create the conversion function. The package will be
  loaded automically if necessary.}

  \item{from}{The ID type to convert from. This should be
  one of the available keytypes in the \code{AnnotationDb}
  package. Use the \code{\link{keytypes}} function to find
  out which keytypes can be used.}

  \item{to}{The ID type to convert to. This should be one
  of the available columns in the \code{AnnotationDb}
  package. Use the \code{\link{cols}} function to find out
  which column names can be used.}
}
\value{
A function which takes a vector of input IDs as single
argument and returns another vector with the converted IDs.
}
\description{
Creates a function based on an \code{\link{AnnotationDb}}
package This package accepts a vector of input IDs and
returns a vector of output IDs. If an input ID cannot be
mapped to an output ID, to output vector will be one
element shorter. This behaviour can be changed by setting
the additional \code{na.rm} argument to \code{FALSE}.
Likewise, if an input ID maps to multiple output IDs, the
output vector will contain all of the latter. If really
need the output vector to have the same length as the input
vector, you can set the \code{drop.ambiguous} argument to
\code{TRUE}
}
\author{
Cedric Simillion
}

