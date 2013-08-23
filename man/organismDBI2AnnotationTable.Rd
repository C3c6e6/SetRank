\name{organismDBI2AnnotationTable}
\alias{organismDBI2AnnotationTable}

\title{Build an annotation frame from the Gene Ontology annotation provided by
organismDBI
}

\description{Generates an annotation data frame that can be used by 
  \code{\link{buildSetCollection}}
}

\usage{
organismDBI2AnnotationTable(annotationPackageName)
}

\arguments{
  \item{annotationPackageName}{The name of the annotationDBI package that
  contains the annotation for the species of interest, e.g. \code{Homo.sapiens}
  or \code{Mus.musculus}.}
}

\value{A data frame containing the annotation data. See 
  \code{\link{buildSetCollection}} for more information.}


\author{Cedric Simillion}

\seealso{
See \code{\link{buildSetCollection}}.
}
\examples{
#Retrieve the GO annotation for human
#annotation = organismDBI2AnnotationTable("Homo.sapiens")
}
\keyword{annotation}
\keyword{gene ontology}
\keyword{GO}
