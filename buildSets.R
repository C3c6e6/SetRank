# Project: SetRank
# 
# Author: cesim
###############################################################################
library("reactome.db")
library("GO.db")
library("igraph")

"%i%" <- intersect
"%u%" <- union
"%d%" <- setdiff

uniqueCount <- function(x) {
	if (class(x) == "factor") length(levels(x)) else length(unique(x))
}

reactome2AnnotationTable <- function(organismName) {
	pathways <- toTable(reactomePATHNAME2ID)
	organismPathways = pathways[grep(organismName,pathways$path_name),]
	organismPathways$path_name = sub(sprintf("^%s: ", organismName), "", 
			organismPathways$path_name)
	rownames(organismPathways) <- organismPathways$reactome_id
	path2GeneID = as.list(reactomePATHID2EXTID)
	nonEmptyPaths = rownames(organismPathways) %i% names(path2GeneID)
	do.call(rbind, lapply(nonEmptyPaths, function(x) data.frame(
								geneID=path2GeneID[[x]], 
								termID=paste("REACTOME", x, sep=":"),
								termName=organismPathways[x,]$path_name, 
								dbName="Reactome", stringsAsFactors=FALSE)))
}

organismDBI2AnnotationTable <- function(annotationPackageName) {
	require(annotationPackageName, character.only=TRUE)
	message("Querying organismDBI...")
	organismDBIData = select(eval(parse(text=annotationPackageName)), 
			keys=keys(eval(parse(text=annotationPackageName)), "GOID"), 
			cols=c("ENTREZID", "TERM"), keytype="GOID")
	message("Constructing preliminary table...")
	organismDBIData = organismDBIData[!is.na(organismDBIData$ENTREZID),]
	preliminaryTable = data.frame(geneID = organismDBIData$ENTREZID, 
			termID = organismDBIData$GOID, termName = organismDBIData$TERM,
			dbName = organismDBIData$ONTOLOGY, stringsAsFactors=FALSE)
	message("Querying GO.db...")
	offspringList = c(as.list(GOBPOFFSPRING), as.list(GOCCOFFSPRING), 
			as.list(GOMFOFFSPRING))
	goTerms = as.list(GOTERM)
	message("Constructing uncovered term table... ", appendLF=FALSE)
	omittedTermIDs = names(offspringList) %d% unique(preliminaryTable$termID)
	message("adding ", length(omittedTermIDs), " omitted terms")
	omittedTermTable = do.call(rbind, lapply(omittedTermIDs, function(x) { 
						t = goTerms[[x]]; 
						data.frame(geneID=NA, termID=GOID(t), termName=Term(t), 
								dbName=Ontology(t), stringsAsFactors=FALSE)}))
	message("merging...")
	preliminaryTable = rbind(preliminaryTable, omittedTermTable)
	message("splitting...")
	tableSplit = split(preliminaryTable, preliminaryTable$termID)
	message("extending...")
	do.call(rbind, lapply(tableSplit, expandWithTermOffspring, tableSplit, 
					offspringList))
}

expandWithTermOffspring <- function(subTable, tableSplit, offspringList) {
	termID = unique(as.character(subTable$termID))
	message(termID)
	termName = unique(as.character(subTable$termName))
	dbName = unique(as.character(subTable$dbName))
	offspring = offspringList[[termID]] %i% names(tableSplit)
	extension = do.call( rbind, lapply(offspring, 
					function(x) data.frame(geneID = tableSplit[[x]]$geneID, 
								termID = termID, termName = termName, 
								dbName = dbName, stringsAsFactors = FALSE)))
	expandedTable = rbind(subTable, extension)
	expandedTable[!is.na(expandedTable$geneID),]
}

buildSetCollection <- function(..., referenceSet = NULL, maxSetSize = 2000) {
	annotationTable = do.call(rbind, list(...))
	if	(!is.null(referenceSet)) {
		annotationTable = 
				annotationTable[annotationTable$geneID %in% referenceSet,]
	} else {
		referenceSet = unique(as.character(annotationTable$geneID))
	}
	collection = list(maxSetSize = maxSetSize, referenceSet=referenceSet)
	collection$sets = by(annotationTable, annotationTable[,"termID"], 
			function(x) {
				geneSet = unique(as.character(x[,"geneID"]))
				attr(geneSet, "ID") <- unique(as.character(x[,"termID"]))
				attr(geneSet, "name") <- unique(as.character(x[,"termName"]))
				attr(geneSet, "db") <- unique(as.character(x[,"dbName"]))
				geneSet
			})
	collection$g =  length(referenceSet)
	setSizes = sapply(collection$sets, length)
	collection$bigSets = names(setSizes[setSizes > maxSetSize])
	message(uniqueCount(annotationTable$dbName), " gene set DBs, ", 
			length(collection$sets), " initial gene sets, ", 
			length(collection$sets) - length(collection$bigSets), 
			" sets remaining and ", collection$g, " genes in collection")
	collection$intersection.p.cutoff = 0.01
	collection$intersections = getSignificantIntersections(collection$sets, 
			annotationTable, collection$g, collection$intersection.p.cutoff,
			collection$bigSets)
	collection
}

getSignificantIntersections <- function(collectionSets, annotationTable, g, 
		pValueCutoff, bigSets) {
	setIDs = names(collectionSets)
	setIndicesPerGene = by(annotationTable, annotationTable$geneID, 
			function(x) which(setIDs %in% 
								(as.character(x$termID) %d% bigSets)))
	setCount = unlist(lapply(setIndicesPerGene, length))
	geneOrder = names(sort(setCount[setCount > 1], decreasing=TRUE))
	intersectionsPerGene = lapply(setIndicesPerGene[geneOrder],
			function(x) apply(t(combn(x, 2)), 1, pack, length(setIDs)))
	intersectionIndices = unlist(intersectionsPerGene, use.names=FALSE)
	intersectionIndices = unique(intersectionIndices)
	message(length(intersectionIndices), " intersections to test...", appendLF=FALSE)
	intersectionPValues = sapply(intersectionIndices, getIntersectionPValue, collectionSets, g)
	significantIndices = which(intersectionPValues <= pValueCutoff)
	pValueFrame = as.data.frame(do.call(rbind, 
					lapply(intersectionIndices[significantIndices], 
					function(x)	setIDs[unpack(x, length(setIDs))])))
	colnames(pValueFrame) <- c("setA", "setB")
	pValueFrame$pValue = intersectionPValues[significantIndices]
	message(nrow(pValueFrame), " intersections significant")
	pValueFrame
}

pack <- function(indexPair, n) {
	if (indexPair[1] > n || indexPair[2] > n) stop("Index higher than n")
	if (indexPair[1] > indexPair[2]) {
		a = indexPair[2]
		b = indexPair[1]
	} else {
		a = indexPair[1]
		b = indexPair[2]
	}
	(a-1)*(n-1)+(b-1)
}

unpack <- function(packed, n) {
	n_ = n-1
	a = ceiling(packed/n_)
	r = (packed %% n_)
	b = if (r == 0) n else r+1
	if (a >= b) stop("Invalid packed value")
	c(a,b)
}

getIntersectionPValue <-function(intersectionIndex, setCollection, g) {
	setIndexPair = unpack(intersectionIndex, length(setCollection))
	setA = setCollection[[setIndexPair[1]]]
	setB = setCollection[[setIndexPair[2]]]
	i = length(setA %i% setB)
	m = length(setA)
	n = length(setB)
	return(intersectionTest(g, m, n, i)$p.value)
}	

intersectionTest <- function(g, m, n, i) {
	fisher.test(rbind(c(i, n-i), c(m, g-(m+n-i))), alternative="greater")
}

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
