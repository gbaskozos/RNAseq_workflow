

GOenrichment <- function(description, geneList, ontology, selection = topDiffGenes, cutoff = 0.05, nodeSize, annType = annFUN.db, annObj, path=getwd(), sig_nodes) {
require(topGO, quietly=TRUE)
require(Rgraphviz, quietly=TRUE)





topDiffGenes <- function (allScore) 
{
    return(allScore < cutoff)
}

cat("Calculating GO values for all genes in ", description, "with a cut off value <", cutoff, ". Using annotation package ", annObj, " and ontology ", ontology, "\n")

godata <- new("topGOdata", description = description, ontology= ontology, allGenes = geneList, geneSel = selection, nodeSize = nodeSize, annotationFun = annType, mapping = annObj, ID = "ensembl")

fisher_weight <- runTest(godata, algorithm = "weight01", statistic = "fisher")

fisher_classic <- runTest(godata, algorithm = "classic", statistic = "fisher")

kselim <- runTest(godata, algorithm = "elim", statistic= "ks")

ksWeight <- runTest(godata, algorithm = "weight01", statistic= "ks")




if(length(score(ksWeight)) >= 10) {
topRes <- GenTable(godata, classicFisher = fisher_classic, weightFisher = fisher_weight, weightKS = ksWeight, elimKS = kselim, orderBy = "weightFisher", ranksOf = "elimKS", topNodes = 10)
write.csv(file = paste(path,"/", description, "_", ontology, "_topGOnodes", ".csv", sep=""), topRes, quote=TRUE, row.names=FALSE)

}

allRes <- GenTable(godata, classicFisher = fisher_classic, weightFisher = fisher_weight, weightKS = ksWeight, elimKS = kselim, orderBy = "weightFisher", ranksOf = "elimKS", topNodes = length(score(fisher_weight)))


write.csv(file = paste(path,"/", description, "_", ontology, "_allGOnodes", ".csv", sep=""), allRes, quote=TRUE, row.names=FALSE)

save(file = paste(path,"/", description, "_", ontology, "_GOdata", ".RData", sep=""), godata)

res_table <- GenTable(godata, classicFisher = fisher_classic, weightFisher = fisher_weight, weightKS = ksWeight, elimKS = kselim, orderBy = "weightFisher", ranksOf = "elimKS", topNodes = length(score(fisher_weight)))

newlist <- list("godata" = godata, "classicFisher" = fisher_classic, "weightFisher" = fisher_weight, "weightKS" = ksWeight, "elimKS" = kselim, "res_sum" = res_table )

return(newlist)
	}

