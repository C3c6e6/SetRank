\name{buildSetCollection}
\alias{buildSetCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a gene set collection}
\description{Builds an object containing the collection of all gene sets to be used 
by the \code{\link{setRankAnalysis}} function.}
\usage{
buildSetCollection(..., referenceSet = NULL, maxSetSize = 2000)
}
\arguments{
  \item{\dots}{One or more data frame objects containing the annotation of 
  genes with pathway identifiers and descriptions. The idea is to provide one
  data frame per pathway database. For Gene Ontology and Reactome annotation,
  you can respectively use the \code{\link{organismDBI2AnnotationTable}} and
  \code{\link{reactome2AnnotationTable}} functions to generate these tables.
  For other pathway databases, you will have to specify your own annotation
  tables. See the \emph{Details} section for more information.
}
  \item{referenceSet}{Optional but recommended. A vector of geneIDs specifying
  the background gene set against which to test for over-representation of gene
  sets. The default is to use all genes present in the supplied gene annotation
  tables. However, many experiments are intrinsically biased for certain
  pathways e.g. because they only contain samples from a specific tissue.
  Supplying a suitable reference set will remove this bias.
}
  \item{maxSetSize}{The maximum number of genes in a gene set. Any gene sets
  with more genes will not be considered during the analysis.
}
}
\details{\subsection{Execution time}{This function typically takes some time to
execute as it pre-calculates all significant intersections between pairs of
gene sets in the collection. An intersection between two gene sets is
considered significant if it contains more elements than expected by chance,
given the sizes of both sets. Performing this calculation beforehand allows to
re-use the same setCollection object for different analysis. It is therefore
recommended to separate the creation of the setCollection object and the actual
analysis in different scripts. Once the collection is created, it can be stored
on disk using the \code{save} command. The analysis script can then load the
collection using the \code{load} command.}
\subsection{Creation of custom annotation tables.}{Any pathway or other gene set
database can be supplied to \code{buildSetCollection}. The only thing that
needs to be done is to provide the annotation data as a data frame containing 
the following fields:
\describe{
  \item{geneID}{The gene identifier. Can be any type of identifier, but one
  must make sure that all annotation frames passed to \code{buildSetCollection}
  use the same identifier. As the output of the 
  \code{\link{organismDBI2AnnotationTable}} and
  \code{\link{reactome2AnnotationTable}} use Entrez Gene identifiers, it is
  best to use these in your own annotation frames as well. Also, make sure the
  identifiers as passed as character and not as integer values.}
  \item{termID}{Pathway identifier. Make sure each pathway identifier is unique
  across all pathway databases used. You can do this by prefixing the IDs with
  a namespace identifier like \code{"REACTOME:"}.}
  \item{termName}{Pathway description. A string describing the pathway, e.g.
  "negative regulation of sterol metabolism"}
  \item{dbName}{A short string given the name of the pathway database used for
  the annotation. E.g. "KEGG".}
}
}
}
\value{A gene set collection which is a list object containing the following
fields:
  \item{maxSetSize}{The maximum set size applied when constructing the
    collection.}
  \item{referenceSet}{A vector listing all gene IDS that are part of the
    reference.}
  \item{sets}{A list of vectors. The list names are the pathway IDs as supplied
    in the \code{termID} column of the annotation frame(s) supplied..}  Each
    vector contains all geneIDs of the gene set and has three attributes set: 
    \code{ID}, \code{name}, and \code{db} which correspond respectively to the
    \code{termID}, \code{termName}, and \code{dbName} fields of the annotation
    frame.
  \item{g}{The size of the reference set.}
  \item{bigSets}{A list of pathway IDs of gene sets with sizes bigger than the
    specified maximum set size.}
  \item{intersection.p.cutoff}{The p-value cutoff used to determine which
    intersections of pairs of gene sets (see \emph{Details}) are significant.}
  \item{intersections}{A data frame listing all significant intersections
    together with the p-value.}  
}
\references{
I'll write a paper about all this! Really! I promise!
}
\author{Cedric Simillion
}

\seealso{
See also \code{\link{organismDBI2AnnotationTable}} and 
\code{\link{reactome2AnnotationTable}}
}

\examples{
#Create a random annotation table and build a set collection from it
randomAnnotation = data.frame(
    geneID = as.character(sample(5000, 20000, replace=TRUE)),
    termID = paste("random", sample(200, 20000, replace=TRUE), sep=":"),
    termName = "random pathway", description="utter randomness",
    dbName = "random")
randomCollection = buildSetCollection(randomAnnotation, maxSetSize=200)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{gene set collection}
\keyword{annotation}
\keyword{set rank analysis}
