\name{setRankAnalysis}
\alias{setRankAnalysis}

\title{Advanced gene set enrichment.
}
\description{Performs advanced gene set enrichment analysis on a set of gene.}
\usage{
setRankAnalysis(setCollection, selectedGenes, setPCutoff = 0.01, heteroPCutoff = setPCutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{setCollection}{A gene set collection object, generated with the
  \code{\link{buildSetCollection}} function.
}
  \item{selectedGenes}{A vector containing the set of gene IDs to test for gene
  set enrichment. This is typically the list of significant genes returned
  by the analysis of an omics dataset.
}
  \item{setPCutoff}{The p-value cutoff to be used to consider a gene set
  significant. Recommended value: 0.01
}
  \item{heteroPCutoff}{The p-value cutoff to used to discover pairs of sets
  with a heterogeneous distribution of significant genes.
}
}

\value{An igraph object. Use the igraph \code{\link{get.data.frame}} function
  to get a data frame with all the significant pathways.
}

\author{
Cedric Simillion
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{buildSetCollection}}
}

\examples{
#Create a random annotation table and build a set collection from it
randomAnnotation = data.frame(
    geneID = as.character(sample(5000, 100000, replace=TRUE)),
    termID = paste("random", sample(200, 100000, replace=TRUE), sep=":"),
    termName = "random pathway",
    dbName = "random")
randomCollection = buildSetCollection(randomAnnotation)

#Create a set of selectedGenes by combining all genes from 10 gene sets with
#1000 random genes and then selecting 2000 from that list.
setSample = sample(names(randomCollection$sets), 10)
selectedGenes = unique(sample(c(unlist(randomCollection$sets[setSample], 
				use.names=FALSE), sample(levels(randomAnnotation$geneID), 
				500)), 2000))
#Perform the gene set enrichment analysis
results = setRankAnalysis(randomCollection, selectedGenes)

}

\keyword{gene set enrichment analysis}
\keyword{set rank}% __ONLY ONE__ keyword per line
