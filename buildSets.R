# Project: SetRank
# 
# Author: cesim
###############################################################################
library("reactome.db")
library("igraph")
uniqueCount <- function(x) {
	if (class(x) == "factor") length(levels(x)) else length(unique(x))
}

"%i%" <- intersect
"%u%" <- union
"%d%" <- setdiff

buildSetCollection <- function(...) {
	annotationTable = do.call(rbind, list(...))
	collection = list()
	collection$sets = by(annotationTable, annotationTable[,"termID"], 
			function(x) {
				geneSet = unique(as.character(x[,"geneID"]))
				attr(geneSet, "ID") <- unique(as.character(x[,"termID"]))
				attr(geneSet, "name") <- unique(as.character(x[,"termName"]))
				attr(geneSet, "db") <- unique(as.character(x[,"dbName"]))
				geneSet
			})
	collection$g = uniqueCount(annotationTable$geneID)
	message(uniqueCount(annotationTable$dbName), " gene set DBs, ", 
			uniqueCount(annotationTable$termID), " gene sets and ", 
			collection$g, " genes in collection")
	collection$intersection.p.cutoff = 0.01
	collection$intersections = getSignificantIntersections(collection$sets, 
			annotationTable, collection$g, collection$intersection.p.cutoff)
	collection
}

organismDBI2AnnotationTable <- function(annotationPackageName) {
	require(annotationPackageName, character.only=TRUE)
	organismDBIData = select(eval(parse(text=annotationPackageName)), 
			keys=keys(eval(parse(text=annotationPackageName)), "GOID"), 
			cols=c("ENTREZID", "TERM"), keytype="GOID")
	organismDBIData = organismDBIData[!is.na(organismDBIData$ENTREZID),]
	data.frame(geneID = organismDBIData$ENTREZID, 
			termID = organismDBIData$GOID, termName = organismDBIData$TERM,
			dbName = organismDBIData$ONTOLOGY)
}

getSignificantIntersections <- function(collectionSets, annotationTable, g, pValueCutoff) {
	setsPerGene = by(annotationTable, annotationTable$geneID, 
			function(x) as.character(x$termID))
	setCount = unlist(lapply(setsPerGene, length))
	geneOrder = names(sort(setCount[setCount > 1], decreasing=TRUE))
	pValueFrame = do.call(rbind, lapply(setsPerGene[geneOrder], 
					function(x) t(combn(x, 2)) ))
	pValueFrame = as.data.frame(unique(pValueFrame, MARGIN=1))
	colnames(pValueFrame) = c("setA", "setB")
	message(nrow(pValueFrame), " intersections to test...", appendLF=FALSE)
	pValueFrame$pValue = apply(pValueFrame, 1, getIntersectionPValue, collectionSets, g)
	pValueFrame = pValueFrame[pValueFrame$pValue <= pValueCutoff,]
	message(nrow(pValueFrame), " intersections significant")
	pValueFrame
}

getIntersectionPValue <-function(setIDPair, setCollection, g) {
	setA = setCollection[[setIDPair[1]]]
	setB = setCollection[[setIDPair[2]]]
	i = length(setA %i% setB)
	m = length(setA)
	n = length(setB)
	return(intersectionTest(g, m, n, i)$p.value)
}	

intersectionTest <- function(g, m, n, i) {
	fisher.test(rbind(c(i, n-i), c(m, g-(m+n-i))), alternative="greater")
}

getPrimarySetPValues <- function(setCollection, selectedGenes, referenceGenes = NA) {
	selectedGenes = as.character(selectedGenes)
	if (!is.na(referenceGenes)) {
		referenceGenes = as.set(as.character(referenceGenes))
		setCollection$sets = lapply(setCollection$sets, intersect, 
				referenceGenes)
		g = length(referenceGenes)
	} else {
		g = setCollection$g
	}
	unlist(lapply(setCollection$sets, getSetPValue, selectedGenes, g))
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
	edgeTable$discardFrom = FALSE
	edgeTable$discardTo = FALSE
	#TODO: check for intersections
	edgeTable[edgeTable$from_pHetero <= heteroPCutoff & 
					edgeTable$from_pDiff > setPCutoff,]$discardFrom = TRUE
	edgeTable[edgeTable$to_pHetero <= heteroPCutoff & 
					edgeTable$to_pDiff > setPCutoff,]$discardTo = TRUE
	edgeTable[edgeTable$discardFrom & 
					edgeTable$discardTrue,]$type = "intersection"
    edgeTable[edgeTable$type == "intersection"]$discardFrom = FALSE
	edgeTable[edgeTable$type == "intersection"]$discardTrue = FALSE
	edgeTable
}

getSetPValue <- function(geneSet, selectedGenes, g) {
	m = length(geneSet)
	i = length(geneSet %i% selectedGenes)
	s = length(selectedGenes)
	fisher.test(rbind(c(i,m-i),c(s-i,g-(m+s-i))), 
			alternative="greater")$p.value
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
	a$pDiff = getSetPValue(diffA, selectedGenes, setCollection$g)
	a$pHetero = setHeterogeneityPValue(diffA, intersection, selectedGenes)
	a$diffSize = length(diffA)
	b$pDiff = getSetPValue(diffB, selectedGenes, setCollection$g)
	b$pHetero = setHeterogeneityPValue(diffB, intersection, selectedGenes)
	b$diffSize = length(diffB)
	if (length(diffA) == 0 || length(diffB) == 0) {
		type = "subset"
		if (length(diffA) == 0) {
			from = a
			to = b
		} else {
			from = b
			to = a
		}
		from$pDiff = 1
		from$pHetero = 1
	} else if (a$pDiff > b$pDiff) {
		from = a
		to = b
	} else {
		from = b
		to = a
	}
	jaccard = length(intersection)/length(unionSet)
	data.frame(from=from$id, to=to$id,  type=type, from_pDiff = from$pDiff, 
			from_pHetero = from$pHetero, from_diffSize = from$diffSize, 
			to_pDiff = to$pDiff, to_pHetero = to$pHetero, 
			to_diffSize = to$diffSize, intersectionSize = length(intersection),
			jaccard=jaccard)
}

setHeterogeneityPValue <- function(difference, intersection, selectedGenes) {
	differenceNotSelection = difference %d% selectedGenes
	differenceSelection = difference %i% selectedGenes
	intersectionNotSelection = intersection %d% selectedGenes
	intersectionSelection = intersection %i% selectedGenes
	fisher.test(rbind(
					c(length(differenceNotSelection),length(differenceSelection)),
					c(length(intersectionNotSelection),length(intersectionSelection))
			))$p.value
}

setRankAnalysis <- function(setCollection, selectedGenes, setPCutoff = 0.01, 
		heteroPCutoff = setPCutoff) {
	pValues = getPrimarySetPValues(setCollection, selectedGenes)
	edgeTable = buildEdgeTable(setCollection, pValues, selectedGenes, 
			setPCutoff, heteroPCutoff)
	toDelete = unique(
			union(edgeTable[discardFrom,]$from, edgeTable[discardTo,]$to))
	vertexTable = data.frame(ID=sapply(setCollection$sets, attr, "ID"), 
			name = sapply(setCollection$sets, attr, "name"), 
			database=sapply(setCollection$sets, attr, "db"),  pValue = pValues,
			size = sapply(setCollection$sets, length))
	vertexTable = vertexTable[vertexTable$pValue <= setPCutoff,]
	message("discarded ", length(toDelete), " out of ", nrow(vertexTable), 
			" gene sets.")
	setNet = graph.data.frame(edgeTable, directed=TRUE, vertices=vertexTable)
	setNet - toDelete
}

load("ratCollection.Rda")
conversionTable = read.table("id_conversion.txt", sep="\t", header=TRUE)
geneIDTable = conversionTable[!is.na(conversionTable$EntrezGene.ID),]
geneAccs = readLines("gene_acc.lst")
geneIDs = unique(c(geneIDTable[geneIDTable$RefSeq.mRNA..e.g..NM_001195597. %in% 
								geneAccs,]$EntrezGene.ID, 
				geneIDTable[geneIDTable$Ensembl.Transcript.ID %in% 
								geneAccs,]$EntrezGene.ID))
testOutput = setRankAnalysis(ratCollection, geneIDs)