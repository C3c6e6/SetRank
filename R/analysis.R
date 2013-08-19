# Project: SetRank
# 
# Author: cesim
###############################################################################
library("igraph")

setRankAnalysis <- function(setCollection, selectedGenes, setPCutoff = 0.01, 
		heteroPCutoff = setPCutoff) {
	selectedGenes = selectedGenes %i% setCollection$referenceSet
	pValues = getPrimarySetPValues(setCollection, selectedGenes)
	edgeTable = buildEdgeTable(setCollection, pValues, selectedGenes, 
			setPCutoff, heteroPCutoff)
	toDelete = unique(union(edgeTable[edgeTable$discardSource,]$source, 
					edgeTable[edgeTable$discardSink,]$sink))
	vertexTable = data.frame(name=sapply(setCollection$sets, attr, "ID"), 
			description = sapply(setCollection$sets, attr, "name"), 
			database=sapply(setCollection$sets, attr, "db"),  pValue = pValues,
			pp = -log10(pValues), size = sapply(setCollection$sets, length))
	vertexTable = vertexTable[vertexTable$pValue <= setPCutoff,]
	message("discarded ", length(toDelete), " out of ", nrow(vertexTable), 
			" gene sets.")
	setNet = graph.data.frame(edgeTable, directed=TRUE, vertices=vertexTable)
	setNet = setNet - toDelete
	subsetEdges = which(E(setNet)$type == "subset")
	setRank = page.rank(setNet - E(setNet)[subsetEdges])
	set.vertex.attribute(setNet, "setRank", index=names(setRank$vector), 
			value=setRank$vector)
}

getPrimarySetPValues <- function(setCollection, selectedGenes) {
	selectedGenes = as.character(selectedGenes)
	g = setCollection$g
	pValues = sapply(setCollection$sets, getSetPValue, selectedGenes, g)
	pValues[setCollection$bigSets] = 1
	pValues
}

getSetPValue <- function(geneSet, selectedGenes, g) {
	m = length(geneSet)
	i = length(geneSet %i% selectedGenes)
	s = length(selectedGenes)
	fisher.test(rbind(c(i,m-i),c(s-i,g-(m+s-i))), 
			alternative="greater")$p.value
}

buildEdgeTable <- function(setCollection, setPValues, selectedGenes, 
		setPCutoff, heteroPCutoff) {
	significantSetIDs = names(setPValues[setPValues <= setPCutoff])
	message(length(significantSetIDs), " significant sets", appendLF=FALSE)
	intersectionTable = setCollection$intersections
	intersectionsToTest =  apply(intersectionTable, 1, 
			function(x) (length(x %i% significantSetIDs) == 2))
	intersectionTable = intersectionTable[intersectionsToTest,]
	message(" - ", nrow(intersectionTable), " intersections to test.")
	edgeTable = do.call(rbind, apply(intersectionTable, 1, getSetPairStatistics, 
					selectedGenes, setCollection))
	edgeTable$discardSource = FALSE
	edgeTable$discardSink = FALSE
	edgeTable[edgeTable$source_pHetero <= heteroPCutoff & 
					edgeTable$source_pDiff > setPCutoff,]$discardSource = TRUE
	edgeTable[edgeTable$sink_pHetero <= heteroPCutoff & 
					edgeTable$sink_pDiff > setPCutoff,]$discardSink = TRUE
	if (any(edgeTable$discardSource & edgeTable$discardSink)) {
		edgeTable[edgeTable$discardSource & 
						edgeTable$discardSink,]$type = "intersection"
		edgeTable[edgeTable$type == "intersection",]$discardSource = FALSE
		edgeTable[edgeTable$type == "intersection",]$discardSink = FALSE
		
	}
	edgeTable
}


getSetPairStatistics <- function(row, selectedGenes, setCollection) {
	setIDA = row[1]
	setIDB = row[2]
	setA = setCollection$sets[[setIDA]]
	setB = setCollection$sets[[setIDB]]
	intersection = setA %i% setB
	a = list(id = setIDA)
	b = list(id = setIDB)
	unionSet = setA %u% setB
	diffA = setA %d% setB
	diffB = setB %d% setA
	type = "overlap"
	pIntersection = getSetPValue(intersection, selectedGenes, setCollection$g)
	a$pDiff = getSetPValue(diffA, selectedGenes, setCollection$g)
	a$pHetero = setHeterogeneityPValue(diffA, setB, selectedGenes)
	a$diffSize = length(diffA)
	a$size = length(setA)
	b$pDiff = getSetPValue(diffB, selectedGenes, setCollection$g)
	b$pHetero = setHeterogeneityPValue(diffB, setA, selectedGenes)
	b$diffSize = length(diffB)
	b$size = length(setB)
	if (length(diffA) == 0 || length(diffB) == 0) {
		type = "subset"
		if (length(diffA) == 0) {
			source = a
			sink = b
		} else {
			source = b
			sink = a
		}
		source$pDiff = 1
		source$pHetero = 1
	} else if (a$pDiff > b$pDiff) {
		source = a
		sink = b
	} else {
		source = b
		sink = a
	}
	jaccard = length(intersection)/length(unionSet)
	data.frame(source=source$id, sink=sink$id, type=type, 
			source_pDiff = source$pDiff, source_ppDiff = -log10(source$pDiff), 
			source_pHetero = source$pHetero, 
			source_ppHetero = -log10(source$pHetero), 
			source_diffSize = source$diffSize, sink_pDiff = sink$pDiff, 
			sink_ppDiff = -log10(sink$pDiff), sink_pHetero = sink$pHetero,
			sink_ppHetero = -log10(sink$pHetero), 
			sink_diffSize = sink$diffSize, 
			deltaP = -log10(sink$pDiff) + log10(source$pDiff),
			intersectionSize = length(intersection), 
			pIntersection = pIntersection, 
			ppIntersection = -log10(pIntersection), jaccard=jaccard, 
			intersectionSourceFraction = length(intersection)/source$size,  
			stringsAsFactors=FALSE)
}

setHeterogeneityPValue <- function(difference, otherSet, selectedGenes) {
	differenceNotSelection = difference %d% selectedGenes
	differenceSelection = difference %i% selectedGenes
	otherSetNotSelection = otherSet %d% selectedGenes
	otherSetSelection = otherSet %i% selectedGenes
	fisher.test(rbind(
					c(length(differenceNotSelection),length(differenceSelection)),
					c(length(otherSetNotSelection),length(otherSetSelection))
			))$p.value
}

