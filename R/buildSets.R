# Project: SetRank
# 
# Author: cesim
###############################################################################

uniqueCount <- function(x) {
	if (class(x) == "factor") length(levels(x)) else length(unique(x))
}

reactome2AnnotationTable <- function(organismName) {
	pathways <- toTable(reactomePATHNAME2ID)
	searchPattern = sprintf("^\\s*%s\\s*:\\s*", organismName)
	organismPathways = unique(pathways[grep(searchPattern,pathways$path_name),])
	organismPathways$path_name = sub(searchPattern, "", 
			organismPathways$path_name)
	rownames(organismPathways) <- organismPathways$DB_ID
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

buildSetCollection <- function(..., referenceSet = NULL, maxSetSize = 500) {
	annotationTable = do.call(rbind, list(...))
	if	(!is.null(referenceSet)) {
		annotationTable = 
				annotationTable[annotationTable$geneID %in% referenceSet,]
		annotationTable$termID = factor(annotationTable$termID)
	} else {
		referenceSet = unique(referenceSet)
		referenceSet = unique(as.character(annotationTable$geneID))
	}
	collection = list(maxSetSize = maxSetSize, referenceSet=referenceSet)
	collection$sets = if (nrow(annotationTable) == 0) list() else
				by(annotationTable, annotationTable[,"termID"], createSet, 
						simplify=FALSE)
	collection$sets[sapply(collection$sets, is.null)] = NULL
	collection$g =  length(referenceSet)
	collection$bigSets = sapply(collection$sets, length) > maxSetSize
	message(uniqueCount(annotationTable$dbName), " gene set DBs, ", 
			length(collection$sets), " initial gene sets, ", 
			length(collection$sets) - length(which(collection$bigSets)), 
			" sets remaining and ", collection$g, " genes in collection")
	collection$intersection.p.cutoff = 0.01
	collection$intersections = getSignificantIntersections(collection$sets, 
			annotationTable, collection$g, collection$intersection.p.cutoff,
			collection$bigSets)
	message("Pre-calculating critical Fisher-test values...")
	collection$iMatrix = fisherCriticalValues(collection$g, maxSetSize, 0.05)
	collection
}

getSignificantIntersections <- function(collectionSets, annotationTable, g, 
		pValueCutoff, bigSets) {
	setIDs = names(collectionSets[!bigSets])
	setIndicesPerGene = by(annotationTable, annotationTable$geneID, 
			function(x) which(setIDs %in% as.character(x$termID)))
	setCount = unlist(lapply(setIndicesPerGene, length))
	geneOrder = names(sort(setCount[(setCount > 1)], decreasing=TRUE))
	intersectionsPerGene = lapply(setIndicesPerGene[geneOrder],
			function(x) apply(t(combn(x, 2)), 1, pack, length(setIDs)))
	intersectionIndices = unlist(intersectionsPerGene, use.names=FALSE)
	intersectionIndices = unique(intersectionIndices)
	message(length(intersectionIndices), " intersections to test...", 
			appendLF=FALSE)
	intersectionPValues = unlist(mclapply(intersectionIndices, 
					getIntersectionPValue, collectionSets[!bigSets], g))
	significantIndices = which(intersectionPValues <= pValueCutoff)
	pValueFrame = as.data.frame(do.call(rbind, 
					lapply(intersectionIndices[significantIndices], 
							function(x)	setIDs[unpack(x, length(setIDs))])),
			stringsAsFactors=FALSE)
	if (nrow(pValueFrame) > 0) {
		colnames(pValueFrame) <- c("setA", "setB")
		pValueFrame$pValue = intersectionPValues[significantIndices]
		message(nrow(pValueFrame), " intersections significant")
	} else {
		pValueFrame = data.frame(setA=c(),setB=c())
	}
	pValueFrame
}

createSet <- function(x) {
	geneSet = unique(as.character(x[,"geneID"]))
	attr(geneSet, "ID") <- 
			unique(as.character(x[,"termID"]))
	attr(geneSet, "name") <- 
			unique(as.character(x[,"termName"]))
	attr(geneSet, "db") <- 
			unique(as.character(x[,"dbName"]))
	attr(geneSet, "description") <- 
			unique(as.character(x[,"description"]))
	geneSet
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

fisherCriticalValues <- function(g,maxSize,criticalP) {
	do.call(rbind, mclapply(1:g, function(s) {
						sapply(1:maxSize, minimalI, s, g, criticalP)
					}))
}

minimalI <- function(m,s,g, maximalP) {
	iVector=min(s,m):0
	pValues = rev(cumsum(dhyper(s-iVector,g-m,m,s)))
	minimalI = which(pValues < maximalP)[1] - 1
}

createIDConverter <- function(annotationPackageName, from, to) {
	require(annotationPackageName, character.only=TRUE)
	keySet = keys(eval(parse(text=annotationPackageName)), from)
	conversionTable = select(eval(parse(text=annotationPackageName)), 
			keys=keySet, keytype=from, columns=to)
	conversion = as.list(by(conversionTable, as.factor(conversionTable[[from]]), 
					function(t) unique(t[[to]]), simplify=FALSE))
	outputFunction = function(x, na.rm=TRUE, drop.ambiguous=FALSE) {
		knownX = unique(x)
		outputList = as.list(rep(NA, length(knownX)))
		names(outputList) = knownX
		knownX = intersect(knownX, names(conversion))
		if (drop.ambiguous) {
			ambiguous = sapply(conversion, function(x) length(x) > 1)
			knownX = setdiff(knownX, names(conversion[ambiguous]))
		}
		outputList[knownX] = conversion[knownX]
		output = unlist(outputList[x], use.names=FALSE)
		if (na.rm) {
			output = output[!is.na(output)]
		}
		output
	} 
	return(outputFunction)
}


