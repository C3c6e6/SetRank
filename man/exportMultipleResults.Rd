% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tables.R
\name{exportMultipleResults}
\alias{exportMultipleResults}
\title{Export multiple SetRank networks and accompanying tables.}
\usage{
exportMultipleResults(networkList, selectedGenesList, collection,
  IDConverter = NULL, outputPath = "./")
}
\arguments{
\item{networkList}{A named list of SetRank networks.}

\item{selectedGenesList}{A named list with the same names as the 
\code{networkList} argument. Each list should be a vector with the set of 
significant genes as Entrez Gene IDs used to construct the SetRank network
with the same name.}

\item{collection}{The set collection used for the SetRank analysis.}

\item{IDConverter}{Optional. By default, Entrez Gene IDs will be displayed
in the output tables. This argument can be used to convert these into more 
human-friendly gene symbols. When supplied, should be a function that takes 
a vector of Entrez Gene IDs as single argument and returns the values of the
corresponding gene symbols or whatever identifier you wish to have displayed
in the output tables.}

\item{outputPath}{The name of the directory where the results should be 
written. If the last element of the path doesn't exist, a directory will be
created.}
}
\value{
None. Files are written out as a side effect.
}
\description{
Given a list of SetRank networks, writes out the following files for each
network:
\enumerate{
\item SetRank network in GML format called \code{<n>.gml} where \code{<n>}
is the network name.
\item A TAB-delimited file listing the signficant pathways in the network,
called \code{<n>_pathways.txt} with \code{<n>} the network name.
\item A TAB-delimited file listing which significant genes belong to which
pathway, called  \code{<n>_membership.txt} with \code{<n>} again the network
name.
}
The network names will be taken from the names of the input list.
Additionally, two global files will be created as well:
\enumerate{
\item A Cytoscape VizMap visualisation file, called\code{setrank.xml}.
\item A TAB-delimited file listing which pathways are found in which
networks, called "pathways.txt".
}
}

