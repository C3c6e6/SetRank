# Project: SetRank
# 
# Author: cesim
###############################################################################
library("reactome.db")
library(sets)

uniqueCount <- function(x) {
	if (class(x) == "factor") length(levels(x)) else length(unique(x))
}

buildSetCollection <- function(...) {
	annotationTable = do.call(rbind, list(...))
	collection = list()
	collection$sets = by(annotationTable, annotationTable[,"termID"], 
			function(x) {
				geneSet = as.set(as.character(x[,"geneID"]))
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
}

setIntersectionsPerGene <- function(setIDs, setCollection, g) {
	message(length(setIDs))
	pValues = unlist(lapply(combn(setIDs, 2, simplify=FALSE), 
					getIntersectionPValue, setCollection, g))
	pValueFrame = as.data.frame(t(combn(setIDs, 2)))
	colnames(pValueFrame) <- c("setA", "setB")
	pValueFrame$pValue <- pValues
	pValueFrame
}

getIntersectionPValue <-function(setIDPair, setCollection, g) {
	setA = setCollection[[setIDPair[1]]]
	setB = setCollection[[setIDPair[2]]]
	i = length(setA & setB)
	m = length(setA)
	n = length(setB)
	return(intersectionTest(g, m, n, i)$p.value)
}	

intersectionTest <- function(g, m, n, i) {
	fisher.test(rbind(c(i, n-i), c(m, g-(m+n-i))), alternative="greater")
}

getPrimarySetPValues <- function(setCollection, selectedGenes, referenceGenes = NA) {
	selectedGenes = as.set(as.character(selectedGenes))
	if (!is.na(referenceGenes)) {
		referenceGenes = as.set(as.character(referenceGenes))
		setCollection$sets = lapply(setCollection$sets, set_intersection, 
				referenceGenes)
		g = length(referenceGenes)
	} else {
		g = setCollection$g
	}
	unlist(lapply(setCollection$sets, getSetPValue, selectedGenes, g))
}

buildEdgeTable <- function(setCollection, setPValues, selectedGenes, setPCutoff = 0.01) {
	significantSetIDs = as.set(names(setPValues[setPValues <= setPCutoff]))
	message(length(significantSetIDs), " significant sets")
	intersectionTable = setCollection$intersections
	intersectionsToTest =  apply(intersectionTable, 1, 
			function(x) length(set_intersection(x, significantSetIDs)) == 2)
	intersectionTable = intersectionTable[intersectionsToTest,]
	do.call(rbind, apply(intersectionTable, 1, getSetPairStatistics, 
					selectedGenes, setCollection))
}

getSetPValue <- function(geneSet, selectedGenes, g) {
	m = length(geneSet)
	i = length(set_intersection(geneSet, selectedGenes))
	s = length(selectedGenes)
	fisher.test(rbind(c(i,m-i),c(s-i,g-(m+s-i))), 
			alternative="greater")$p.value
}

getSetPairStatistics <- function(row, selectedGenes, setCollection) {
	setIDA = row[1]
	setIDB = row[2]
	message(setIDA, " ", setIDB)
	setA = setCollection$sets[[setIDA]]
	setB = setCollection$sets[[setIDB]]
	intersection = setA & setB
	a = list(id = setIDA)
	b = list(id = setIDB)
	union = setA | setB
	diffA = setA - setB
	diffB = setB - setA
	type = "intersection"
	a$pDiff = getSetPValue(diffA, selectedGenes, setCollection$g)
	a$pHetero = setHeterogeneityPValue(diffA, intersection, selectedGenes)
	b$pDiff = getSetPValue(diffB, selectedGenes, setCollection$g)
	b$pHetero = setHeterogeneityPValue(diffB, intersection, selectedGenes)
	if (length(diffA) == 0 || length(diffB) == 0) {
		type = "subset"
		if (length(diffA) == 0) {
			donor = a
			acceptor = b
		} else {
			donor = b
			acceptor = a
		}
		donor$pDiff = 1
		donor$pHetero = 1
	} else if (a$pDiff > b$pDiff) {
		donor = a
		acceptor = b
	} else {
		donor = b
		acceptor = a
	}
	jaccard = length(intersection)/length(union)
	data.frame(donor=donor$id, type=type, acceptor=acceptor$id, 
			donor_pDiff = donor$pDiff, donor_pHetero = donor$pHetero, 
			acceptor_pDiff = acceptor$pDiff, 
			acceptor_pHetero = acceptor$pHetero, jaccard=jaccard)
}

setHeterogeneityPValue <- function(difference, intersection, selectedGenes) {
	differenceNotSelection = difference - selectedGenes
	differenceSelection = difference & selectedGenes
	intersectionNotSelection = intersection - selectedGenes
	intersectionSelection = intersection & selectedGenes
	fisher.test(rbind(
		c(length(differenceNotSelection),length(differenceSelection)),
		c(length(intersectionNotSelection),length(intersectionSelection))
   ))$p.value
}

load("ratCollection.Rda")
conversionTable = read.table("id_conversion.txt", sep="\t", header=TRUE)
geneIDTable = conversionTable[!is.na(conversionTable$EntrezGene.ID),]
geneAccs = readLines("gene_acc.lst")
geneIDs = unique(c(geneIDTable[geneIDTable$RefSeq.mRNA..e.g..NM_001195597. %in% 
								geneAccs,]$EntrezGene.ID, 
				geneIDTable[geneIDTable$Ensembl.Transcript.ID %in% 
								geneAccs,]$EntrezGene.ID))
message(Sys.time())
pValues = getPrimarySetPValues(ratCollection, geneIDs)
message(Sys.time())
testTable = buildEdgeTable(ratCollection, pValues, geneIDs)
message(Sys.time())
