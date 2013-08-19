# Project: SetRank
# 
# Author: cesim
###############################################################################


library(hash)
library(sets)
library(ggplot2)

hash_intersect <- function(hashA, hashB) {
	commonIndex = has.key(keys(hashB), hashA)
	commonKeys = names(commonIndex[commonIndex])
	return(hashA[commonKeys])
}

hash_union <- function(hashA, hashB) {
	hash(keys=c(keys(hashA), keys(hashB)), values=TRUE)
}

hash_difference <- function(hashA, hashB) {
	commonIndex = has.key(keys(hashA), hashB)
	differenceKeys = names(commonIndex[!commonIndex])
	return(hashA[differenceKeys])
}

list_intersect <- function(listA, listB) {
	intersection = list()
	for (keyA in names(listA)) {
		if (!is.null(listB[[keyA]])) {
			intersection[[keyA]] = TRUE
		}
	}
	intersection
}

list_union <- function(listA, listB) {
	unionList = listB
	for (keyA in names(listA)) {
		if (is.null(listB[[keyA]])) {
			unionList[[keyA]] = TRUE
		}
	}
	unionList
}

list_difference <- function(listA, listB) {
	difference = listA
	for (keyA in names(listA)) {
		if (!is.null(listB[[keyA]])) {
			difference[[keyA]] = NULL
		}
	}
	difference
}


randomSet <- function(size) {
	as.set(as.character(sample(1:(size*3), size)))
}

randomHashSet <- function(size) {
	hash(keys=as.character(sample(1:(size*3), size)), values=TRUE)
}

randomListSet <- function(size) {
	listSet = as.list(rep(TRUE, size))
	names(listSet) <- as.character(sample(1:(size*3), size))
	listSet
}

randomVectorSet <- function(size) {
	as.character(sample(1:(size*3), size))
}

testSets <- function(size, reps) {
	for (i in 1:reps) {
		a = randomSet(size)
		b = randomSet(size)
		intersectionSet = set_intersection(a,b)
		unionSet = set_union(a,b)
		differenceSet = a - b
	}
}

testHash <- function(size, reps) {
	for (i in 1:reps) {
		a = randomHashSet(size)
		b = randomHashSet(size)
		intersectionHash = hash_intersect(a,b)
		unionHash = hash_union(a,b)
		differenceHash = hash_difference(a,b)
	}
}

testList <- function(size, reps) {
	for (i in 1:reps) {
		a = randomListSet(size)
		b = randomListSet(size)
		names(b) <- as.character(sample(1:(size*3), size))
		intersectionList = list_intersect(a,b)
		unionList = list_union(a,b)
		differenceList = list_difference(a,b)	
	}
}

testVector <- function(size, reps) {
	for (i in 1:reps) {
		a = randomVectorSet(size)
		b = randomVectorSet(size)
		intersectionVector = intersect(a,b)
		unionVector = union(a,b)
		differenceVector = setdiff(a,b)
	}
}

timeTest <- function(size, reps, testFunction) {
	message(size)
	system.time(testFunction(size, reps))[1]
}

sizes = 2^(1:14)
replicates = 10

vectorTotalTimes = sapply(sizes, timeTest, reps = replicates, testVector)
listTotalTimes = sapply(sizes, timeTest, reps = replicates, testList)
setTotalTimes = sapply(sizes, timeTest, reps = replicates, testSets)
hashTotalTimes = sapply(sizes, timeTest, reps = replicates, testHash)

setConstructionTimes = sapply(sizes, function(s) 
			system.time(for(i in 1:(replicates*2)) randomSet(s))[1])
hashConstructionTimes = sapply(sizes, function(s) 
			system.time(for(i in 1:(replicates*2)) randomHashSet(s))[1])
listConstructionTimes = sapply(sizes, function(s) 
			system.time(for(i in 1:(replicates*2)) randomListSet(s))[1])
vectorConstructionTimes = sapply(sizes, function(s) 
			system.time(for(i in 1:(replicates*2)) randomVectorSet(s))[1])

listData = data.frame(size=sizes, totalTime = listTotalTimes, 
		constructionTime=listConstructionTimes, method="list")
setData = data.frame(size=sizes, totalTime = setTotalTimes, 
		constructionTime=setConstructionTimes, method="set")
hashData = data.frame(size=sizes, totalTime = hashTotalTimes, 
		constructionTime=hashConstructionTimes, method="hash")
vectorData = data.frame(size=sizes, totalTime = vectorTotalTimes, 
		constructionTime=vectorConstructionTimes, method="vector")

allData = rbind(listData, setData, hashData, vectorData)
allData$executionTime = allData$totalTime - allData$constructionTime

p <- ggplot(data = allData, aes(x = size, y = executionTime))
p <- p + labs(x = "set size")
p <- p + labs(y = "execution time (s)")
p <- p + geom_line(data = allData, 
		aes(x = size, y = executionTime, group = method, colour= method))
p <- p + geom_point(data = allData, 
		aes(x = size, y = executionTime, group = method, colour= method))
print(p)