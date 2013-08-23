\name{reactome2AnnotationTable}
\alias{reactome2AnnotationTable}
\title{
Build an annotation frame from the annotation data provided by the 
\code{reactome.db} package
}
\description{Converts the \code{reactome.db} data into an annotation data frame
  that can be used by \code{\link{buildSetCollection}}
}

\usage{
reactome2AnnotationTable(organismName)
}

\arguments{
  \item{organismName}{A character value listing the name of the organism for
  which to retrieve the Reactome annotation, e.g. \code{"Homo sapiens"} or 
  \code{"Mus musculus"}
}
}
\value{A data frame containing the annotation data. See 
  \code{\link{buildSetCollection}} for more information.
}

\author{Cedric Simillion}

\seealso{
See \code{\link{buildSetCollection}}.
}

\examples{
#Retrieve the Reactome annotation for human:
annotation = reactome2AnnotationTable("Homo sapiens")
}

\keyword{Reactome}
\keyword{annotation}
