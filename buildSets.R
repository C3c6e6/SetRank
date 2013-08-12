# Project: SetRank
# 
# Author: cesim
###############################################################################
library("reactome.db")
library("go.db")

library("igraph")
uniqueCount <- function(x) {
	if (class(x) == "factor") length(levels(x)) else length(unique(x))
}

"%i%" <- intersect
"%u%" <- union
"%d%" <- setdiff

buildSetCollection <- function(..., maxSetSize = 2000) {
	annotationTable = do.call(rbind, list(...))
	collection = list(maxSetSize = maxSetSize)
	collection$sets = by(annotationTable, annotationTable[,"termID"], 
			function(x) {
				geneSet = unique(as.character(x[,"geneID"]))
				attr(geneSet, "ID") <- unique(as.character(x[,"termID"]))
				attr(geneSet, "name") <- unique(as.character(x[,"termName"]))
				attr(geneSet, "db") <- unique(as.character(x[,"dbName"]))
				geneSet
			})
	collection$g = uniqueCount(annotationTable$geneID)
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

getSignificantIntersections <- function(collectionSets, annotationTable, g, 
		pValueCutoff, bigSets) {
	setsPerGene = by(annotationTable, annotationTable$geneID, 
			function(x) as.character(x$termID) %d% bigSets)
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
	pValues = sapply(setCollection$sets, getSetPValue, selectedGenes, g)
	pValues[setCollection$bigSets] = 1
	pValues
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
	edgeTable[edgeTable$discardSource & 
					edgeTable$discardSink,]$type = "intersection"
    edgeTable[edgeTable$type == "intersection",]$discardSource = FALSE
	edgeTable[edgeTable$type == "intersection",]$discardSink = FALSE
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
    pIntersection = getSetPValue(intersection, selectedGenes, setCollection$g)
	a$pDiff = getSetPValue(diffA, selectedGenes, setCollection$g)
	a$pHetero = setHeterogeneityPValue(diffA, intersection, selectedGenes)
	a$diffSize = length(diffA)
	a$size = length(setA)
	b$pDiff = getSetPValue(diffB, selectedGenes, setCollection$g)
	b$pHetero = setHeterogeneityPValue(diffB, intersection, selectedGenes)
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
			sink_diffSize = sink$diffSize, 
			deltaP = -log10(sink$pDiff) + log10(source$pDiff),
			intersectionSize = length(intersection), 
			pIntersection = pIntersection, 
			ppIntersection = -log10(pIntersection), jaccard=jaccard, 
			intersectionSourceFraction = length(intersection)/source$size,  
			stringsAsFactors=FALSE)
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
	setNet - toDelete
}

#ratTable = organismDBI2AnnotationTable("Rattus.norvegicus")
#load("ratTable.Rda")
#ratCollection=buildSetCollection(ratTable, maxSetSize = 500)

load("ratCollection.Rda")
conversionTable = read.table("id_conversion.txt", sep="\t", header=TRUE)
geneIDTable = conversionTable[!is.na(conversionTable$EntrezGene.ID),]
geneAccs = readLines("gene_acc.lst")
geneIDs = unique(c(geneIDTable[geneIDTable$RefSeq.mRNA..e.g..NM_001195597. %in% 
								geneAccs,]$EntrezGene.ID, 
				geneIDTable[geneIDTable$Ensembl.Transcript.ID %in% 
								geneAccs,]$EntrezGene.ID))
testOutput = setRankAnalysis(ratCollection, geneIDs)
write.graph(testOutput, "ratBrain.gml", format="gml")