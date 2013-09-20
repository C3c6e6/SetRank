# Project: SetRank
# 
# Author: cesim
###############################################################################
library("igraph")

setRankAnalysis <- function(setCollection, selectedGenes, setPCutoff = 0.01, 
		fdrCutoff = 0.05) {
	selectedGenes = selectedGenes %i% setCollection$referenceSet
	pValues = getPrimarySetPValues(setCollection, selectedGenes)
	edgeTable = buildEdgeTable(setCollection, pValues, selectedGenes, 
			setPCutoff)
	toDelete = getNodesToDelete(edgeTable)
	vertexTable = data.frame(name=sapply(setCollection$sets, attr, "ID"), 
			description = sapply(setCollection$sets, attr, "name"), 
			database=sapply(setCollection$sets, attr, "db"),  pValue = pValues,
			pp = -log10(pValues), size = sapply(setCollection$sets, length),
			stringsAsFactors=FALSE)
	vertexTable$nSignificant = sapply(setCollection$sets, 
			function(x) length(x %i% selectedGenes))
	vertexTable = vertexTable[vertexTable$pValue <= setPCutoff,]
	message("discarded ", length(toDelete), " out of ", nrow(vertexTable), 
			" gene sets.")
	setNet = if (is.na(edgeTable[1,]$source)) {
				add.vertices(graph.empty(), nrow(vertexTable),
						attr=as.list(vertexTable))
			} else {
				graph.data.frame(edgeTable, directed=TRUE, 
						vertices=vertexTable)
			}
	if (length(toDelete) > 0) setNet = setNet - toDelete
	subsetEdges = which(E(setNet)$type == "subset")
	setRank = page.rank(setNet - E(setNet)[subsetEdges])
	setNet = set.vertex.attribute(setNet, "setRank", index=names(setRank$vector), 
			value=setRank$vector)
	setNet = addAdjustedPValues(setNet)
	notSignificant = which(V(setNet)$adjustedPValue > fdrCutoff)
	message(length(notSignificant), " gene sets removed after FDR correction.")
	setNet - V(setNet)[notSignificant]
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
		setPCutoff) {
	significantSetIDs = names(setPValues[setPValues <= setPCutoff])
	message(length(significantSetIDs), " significant sets", appendLF=FALSE)
	intersectionTable = setCollection$intersections
	intersectionsToTest =  apply(intersectionTable, 1, 
			function(x) (length(x %i% significantSetIDs) == 2))
	intersectionTable = intersectionTable[intersectionsToTest,]
	message(" - ", nrow(intersectionTable), " intersections to test.")
	if (nrow(intersectionTable) == 0) {
		return(data.frame(source=NA, sink=NA, type=NA))
	}
	edgeTable = do.call(rbind, apply(intersectionTable, 1, getSetPairStatistics, 
					selectedGenes, setCollection))
	edgeTable$discardSource = FALSE
	edgeTable$discardSink = FALSE
	discardSourceIndices = edgeTable$source_pDiff > setPCutoff & 
			edgeTable$type != "subset"
	if (any(discardSourceIndices)) {
		edgeTable[discardSourceIndices,]$discardSource = TRUE
	}
	discardSinkIndices = edgeTable$sink_pDiff > setPCutoff
	if (any(discardSinkIndices)) {
		edgeTable[discardSinkIndices,]$discardSink = TRUE
	}
	if (any(edgeTable$discardSource & edgeTable$discardSink)) {
		edgeTable[edgeTable$discardSource & 
						edgeTable$discardSink,]$type = "intersection"
		edgeTable[edgeTable$type == "intersection",]$discardSource = FALSE
		edgeTable[edgeTable$type == "intersection",]$discardSink = FALSE
	}
	edgeTable = edgeTable[edgeTable$significantJaccard > 0,]
	edgeTable
}


getSetPairStatistics <- function(row, selectedGenes, setCollection) {
	setIDA = row[1]
	setIDB = row[2]
	setA = setCollection$sets[[setIDA]]
	setB = setCollection$sets[[setIDB]]
	intersection = setA %i% setB
	intersectionSignificant = length(intersection %i% selectedGenes)
	a = list(id = setIDA, size = length(setA))
	b = list(id = setIDB, size = length(setB))
	unionSet = setA %u% setB
	diffA = setA %d% setB
	diffB = setB %d% setA
	type = "overlap"
	pIntersection = getSetPValue(intersection, selectedGenes, setCollection$g)
	a$pDiff = getSetPValue(diffA, selectedGenes, setCollection$g)
	a$diffSize = length(diffA)
	a$diffSignificant = length(diffA %i% selectedGenes)
	b$pDiff = getSetPValue(diffB, selectedGenes, setCollection$g)
	b$diffSize = length(diffB)
	b$diffSignificant = length(diffB %i% selectedGenes)
	if (length(diffA) == 0 || length(diffB) == 0) {
		type = "subset"
		if (length(diffA) == 0) {
			source = a
			sink = b
		} else {
			source = b
			sink = a
		}
		source$pDiff = 1.0
	} else if (a$pDiff > b$pDiff) {
		source = a
		sink = b
	} else {
		source = b
		sink = a
	}
	jaccard = length(intersection)/length(unionSet)
	significantJaccard = length(intersection %i% selectedGenes) /
			length(unionSet %i% selectedGenes)
	data.frame(source=source$id, sink=sink$id, type=type, 
			source_pDiff = source$pDiff, source_ppDiff = -log10(source$pDiff), 
			source_diffSize = source$diffSize, 
			source_diffSignificant = source$diffSignificant, 
			sink_pDiff = sink$pDiff, sink_ppDiff = -log10(sink$pDiff), 
			sink_diffSize = sink$diffSize, 
			sink_diffSignificant = sink$diffSignificant,
			deltaP = -log10(sink$pDiff) + log10(source$pDiff),
			intersectionSize = length(intersection), 
			intersectionSignificant = intersectionSignificant,
			pIntersection = pIntersection, 
			ppIntersection = -log10(pIntersection), jaccard=jaccard, 
			significantJaccard = significantJaccard,
			intersectionSourceFraction = length(intersection)/source$size,  
			stringsAsFactors=FALSE)
}

getNodesToDelete <- function(edgeTable) {
	if (is.na(edgeTable[1,]$sink)) {return(c())}
	superSetsToDelete = edgeTable[edgeTable$type == "subset" & 
					edgeTable$discardSink,]$sink
	killTable = edgeTable[edgeTable$type == "overlap",]
	while (TRUE) {
		killTable = killTable[killTable$discardSource,]
		allKillers = unique(killTable$sink)
		allKilled = unique(killTable$source)
		topKillers = allKillers %d% allKilled
		intermediateKillers = 
				unique(killTable[killTable$sink %in% topKillers,]$source) %i%
				allKillers
		invalidKills = killTable$sink %in% intermediateKillers
		if (any(invalidKills)) {
			killTable[invalidKills,]$discardSource = FALSE
		} else {
			break;
		}
	}
	return(unique(superSetsToDelete %u% killTable$source))
}

fisherTest <- function(difference, otherSet, selectedGenes) {
	differenceNotSelection = difference %d% selectedGenes
	differenceSelection = difference %i% selectedGenes
	otherSetNotSelection = otherSet %d% selectedGenes
	otherSetSelection = otherSet %i% selectedGenes
	fisher.test(rbind(
					c(length(differenceNotSelection),length(differenceSelection)),
					c(length(otherSetNotSelection),length(otherSetSelection))
			))$p.value
}

binomialIntervalTest <- function(difference, otherSet, selectedGenes) {
	if (length(difference) == 0) {
		return(1)
	}
	differenceSelection = difference %i% selectedGenes
	otherSetSelection = otherSet %i% selectedGenes
	diffInterval = binomialInterval(
			length(differenceSelection)/length(difference), 
			length(difference), 0.01)
	otherInterval = binomialInterval(length(otherSetSelection)/length(otherSet),
			length(otherSet), 0.01)
	if ((diffInterval$upper > otherInterval$upper) && 
			(diffInterval$lower <= otherInterval$upper)) {
		return(1)
	} else if ((otherInterval$upper > diffInterval$upper) && 
			(otherInterval$lower <= diffInterval$upper)) {
		return(1)
	} else {
		return(0)
	}
}

binomialInterval <- function(p, n, alfa) {
	z = qnorm(1-alfa/2)
	term1 = 2*n*p + z^2
	term2 = z * sqrt( z^2 - 1/n + 4*n*p*(1-p) + (4*p-2) ) + 1
	denominator = 2 * (n + z^2)
	lower = max(0, (term1-term2)/denominator)
	term2 = z * sqrt( z^2 - 1/n + 4*n*p*(1-p) - (4*p-2) ) + 1
	upper = min(1, (term1+term2)/denominator)
	return(list(lower=lower, upper=upper))
}

setHeterogeneityPValue <- fisherTest

addAdjustedPValues <- function(setNet) {
	if (vcount(setNet) == 0) {
		return(setNet)
	}
	edgeTable = get.data.frame(setNet, what="edges")
	nodeTable = get.data.frame(setNet, what="vertices")
	nodeTable$correctedPValue = nodeTable$pValue
	if (ecount(setNet) > 0) {
		correctedPValues = getCorrectedPValues(edgeTable)
		nodeTable[names(correctedPValues),]$correctedPValue = correctedPValues
	}
	nodeTable$adjustedPValue = p.adjust(nodeTable$correctedPValue)
	nodeTable$pp = -log10(nodeTable$correctedPValue)
	setNet = set.vertex.attribute(setNet, "correctedPValue", 
			index=rownames(nodeTable), value=nodeTable$correctedPValue)
	setNet = set.vertex.attribute(setNet, "adjustedPValue", 
			index=rownames(nodeTable), value=nodeTable$adjustedPValue)
	setNet = set.vertex.attribute(setNet, "pp", index=rownames(nodeTable),
			value=nodeTable$pp)
	return(setNet)
}

getMaxP <- function(edgeTable, type, idAttribute, pAttribute) {
	subTable = edgeTable[edgeTable$type == type,]
	maxPList = by(subTable, as.factor(subTable[[idAttribute]]),
			function(x) max(x[[pAttribute]]), simplify=FALSE)
	unlist(maxPList)
}

getCorrectedPValues <- function(edgeTable) {
	maxPValues = list(
			subset = getMaxP(edgeTable, "subset", "to", "sink_pDiff"),
			overlap = getMaxP(edgeTable, "overlap", "from", "source_pDiff"),
			intersectFrom = getMaxP(edgeTable, "intersection", "from", 
					"pIntersection"),
			intersectTo = getMaxP(edgeTable, "intersection", "to", 
					"pIntersection"))
	setIDs = unique(unlist(lapply(maxPValues, names), use.names=FALSE))
	maxPTable = do.call(cbind, lapply(maxPValues, function(x) x[setIDs]))
	rownames(maxPTable) <- setIDs
	return(apply(maxPTable, 1, max, na.rm=TRUE))
}
