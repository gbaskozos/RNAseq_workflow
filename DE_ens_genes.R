library("GenomicRanges")
library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library(genefilter)
library("vsn")
library("BiocParallel")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("sva")
library(Tmisc)
library(calibrate)
library("biomaRt")
library("IHW")
library(lattice)
library(plotrix)
library(WGCNA)

n_cores <- detectCores() - 1

multicoreParam <- MulticoreParam(workers = n_cores, progressbar = TRUE)
options(mc.cores=n_cores)
register(multicoreParam)

allowWGCNAThreads(n_cores)

options(stringsAsFactors=FALSE)

registered()

PATH = "/home/george/Desktop/arc_ndcn_data/RNAseq/star_out"

system("mkdir /home/george/Desktop/arc_ndcn_data/RNAseq/DE") 

PATH_results = "/home/george/Desktop/arc_ndcn_data/RNAseq/DE"

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

########################################
##Create file lists for each condition##
########################################
count_files <- list.files(path="/home/george/Desktop/arc_ndcn_data/RNAseq/star_out", pattern="*.bam.counts.tab$", full.names=TRUE)

file_names <- list.files(path="/home/george/Desktop/arc_ndcn_data/RNAseq/star_out", pattern="*.bam.counts.tab$", full.names=TRUE)

bam_names <- substr(file_names, 69, nchar(file_names)-40)

colData <- data.frame(condition = c(rep("WT",4), rep("KO",4)), id = bam_names, names = c("AD1", "AD2", "AD3", "AD4", "K1", "K2", "K3", "K4"))


###########################
colData$condition <- factor(colData$condition)


toc <- lapply(count_files,read.table, sep="\t", row.names=1)


toc <- data.frame(toc[1:length(bam_names)]) 

names(toc) <- paste(colData$id, colData$names, sep="_")

no_counts <- tail(toc)[-1,]
toc <- toc[1:(dim(toc)[1]-5),]

write.csv(toc, paste0(PATH_results, "/ens_genes_toc_all.csv"))

save(toc, file = paste0(PATH_results, "/ens_genes_toc_all.RData"))

###################
# Calculate TPMs #
###################
library(biomaRt)
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
transcript_lengths <- getBM(attributes=c('ensembl_gene_id', 'transcript_length'), filters = 'ensembl_gene_id', values = rownames(toc), mart = human)
gene_length <- aggregate(transcript_length ~ ensembl_gene_id, data = transcript_lengths, FUN = mean)

toc_for_tpm <- toc[rownames(toc) %in% gene_length$ensembl_gene_id,]

toc_TPM <- apply(toc_for_tpm, 2, tpm, lengths=gene_length$transcript_length)

cut_off <- 0
toc_TPM <- toc_TPM[apply(toc_TPM, 1, function(x) sum(x[1:ncol(toc_TPM)] > cut_off) >= 1), ]

#write.csv(file = paste0(PATH_results, "/TPM_wt_vs_ko.csv"), toc_TPM)
save(file = paste0(PATH_results, "/TPM_wt_vs_ko.csv"),toc_TPM)
##########
# DESEq2 #
##########
dds <- DESeqDataSetFromMatrix(countData = toc, colData = colData, design = ~ condition)
dds$condition <- relevel(dds$condition,"WT")
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
norm_counts <- counts(dds, normalized=TRUE)
colnames(norm_counts) <- names(toc)

pdf(file = paste0(PATH_results,"/SCN9A_counts.pdf"))
plotCounts(dds, gene="ENSG00000169432", intgroup="condition", main = "SCN9A counts")
dev.off()

###################
# Calculate FPKMs #
###################
rownames(gene_length) <- gene_length$ensembl_gene_id

gene_length <- gene_length[rownames(gene_length) %in% rownames(dds),]

dds_fpkm <- dds[rownames(dds) %in% rownames(gene_length)]

mcols(dds_fpkm)$basepairs <- gene_length$transcript_length

KO_WT_fpkm <- fpkm(dds_fpkm, robust=TRUE)
save(file = paste0(PATH_results, "/KO_vs_wt_fpkm.RData"), KO_WT_fpkm)


####################
# Sample distances #
####################
rld_mat <- assay(rld)
colnames(rld_mat) <- names(toc)

save(rld_mat, rld, norm_counts, file = paste(PATH_results, "/ens_genes_norm_counts.RData", sep=""))

hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)


distsRL <- dist(t(rld_mat))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- names(toc)
hc <- hclust(distsRL)

pdf(file = paste(PATH_results, "/heatmap_samples.pdf", sep=""))
heatmap.2(mat, col = rev(hmcol), Rowv=as.dendrogram(hc), trace="none", symm=TRUE, margin=c(10,10), cexCol=0.4, cexRow=0.4)
dev.off()


pdf(file = paste(PATH_results, "/pca_samples.pdf", sep=""))
d = plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(d, "percentVar"))
ggplot(d, aes(x=PC1,y=PC2, color=condition, label=names(toc))) + geom_point(size=4) + geom_text(size=1, check_overlap=TRUE, vjust = 0, nudge_y = 0.5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()

##########################

  rv <- rowVars(assay(vsd))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(vsd)[select,]))


pca_1 <- pca$rotation[,1]
pca_1 <- pca_1[order(abs(pca_1))]
pca_1_symbol <- mapIds(org.Hs.eg.db, keys=names(pca_1), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

pca_2 <- pca$rotation[,2]
pca_2 <- pca_2[order(abs(pca_2))]
pca_2_symbol <- mapIds(org.Hs.eg.db, keys=names(pca_2), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

top_pca <- data.frame(PC1=names(pca_1)[1:50], PC1_symbol = pca_1_symbol[1:50], PC2=names(pca_2)[1:50], PC2=pca_2_symbol[1:50], row.names=NULL)

write.csv(file = paste(PATH_results, "/top_pca.csv", sep=""), top_pca)



#########################################
save(rld, vsd, dds, colData, bam_names, toc, count_files, file = paste(PATH_results, "/DE_genes.RData", sep=""))

###############
# DE results  #
###############
########
#DESeq2#
########

####with Beta prior
dds <- nbinomWaldTest(dds, betaPrior=TRUE)
resultsNames(dds)
#unweighted
res <- results(dds, contrast=c("condition","KO","WT"))
save(file = paste(PATH_results, "/DDS_betaPrior.RData", sep=""), dds, res)

res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
write.csv(file = paste(PATH_results, "/res_all.csv" , sep=""), res[order(-abs(res$log2FoldChange)),])


sig_res <- res[res$padj < 0.05 & !is.na(res$padj),]

write.csv(file = paste(PATH_results, "/sig_de.csv" , sep=""), sig_res[order(-abs(sig_res$log2FoldChange), sig_res$padj),])

sig_res <- sig_res[order(-abs(sig_res$log2FoldChange), sig_res$padj),]

#################
# Create plots #
#################

pdf(file = paste(PATH_results, "/volcano_plot.pdf", sep=""))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, col="grey", main="Volcano plot KO vs WT"))
legend("topleft", c("Adjusted Pvalue < 0.05","Log2 Fold Change > 1", "Adj Pvalue < 0.05 & Lfc > 1"), pch=c(20,20,20), col=c("green","orange", "red"))
with(subset(res, padj< .05),points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(res, abs(log2FoldChange)>1),points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1),points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(sig_res[c(1:50),], textxy(log2FoldChange, -log10(pvalue), labs=symbol, cex=.3))
dev.off()

rld_eset <- as.data.frame(rld_mat)
rld_eset$symbol <- mapIds(org.Hs.eg.db, keys=row.names(rld_eset), column="SYMBOL", keytype="ENSEMBL", multiVals="first")


source("GOenrichment.R")
source("plot_go.R")

cut_off <- 0
toc_cutoff <- toc[apply(toc, 1, function(x) sum(x[1:ncol(toc)] > cut_off) == ncol(toc)), ]

res_list <- res[rownames(toc_cutoff), ]$padj
names(res_list) <- rownames(res[rownames(toc_cutoff), ])
res_list <- res_list[!is.na(res_list)]
res_list_BP <- GOenrichment("KO_vs_WT", geneList = res_list, ontology = "BP", cutoff = 0.05, nodeSize = 10, annType = annFUN.org, annObj = "org.Hs.eg.db", path = PATH_results, sig_nodes=5)
pdf(file = paste(PATH_results, "KO_vs_WT", "_", "BP", "_sigGO_graph",".pdf", sep=""))
showSigOfNodes(GOdata=res_list_BP$godata, termsP.value = score(res_list_BP$weightFisher), plotFunction = GOplot_custom, firstSigNodes = 5, useInfo = 'all', .NO.CHAR = 47, sigForAll=TRUE, sigFontSize = 70, sigBoxSize = c(36,31,31,14,14), elseFontSize = 18)
dev.off()

select_Genes <- function (allScore) 
{
    return(allScore == 1)
}

#Upregulated
res_up <- res[res$log2FoldChange > 0 & !is.na(res$log2FoldChange),]
res_list <- res_up[rownames(toc_cutoff), ]$padj
names(res_list) <- rownames(res_up[rownames(toc_cutoff), ])
res_list <- res_list[!is.na(res_list)]
res_list_BP <- GOenrichment("KO_vs_WT_upregulated", geneList = res_list, ontology = "BP", cutoff = 0.05, nodeSize = 10, annType = annFUN.org, annObj = "org.Hs.eg.db", path = PATH_results, sig_nodes=5)

pdf(file = paste(PATH_results, "/KO_vs_WT_upregulated", "_", "BP", "_sigGO_graph",".pdf", sep=""))
showSigOfNodes(GOdata=res_list_BP$godata, termsP.value = score(res_list_BP$weightFisher), plotFunction = GOplot_custom, firstSigNodes = 5, useInfo = 'all', .NO.CHAR = 47, sigForAll=TRUE, sigFontSize = 70, sigBoxSize = c(36,31,31,14,14), elseFontSize = 18)
dev.off()

#Downregulated
res_down <- res[res$log2FoldChange < 0 & !is.na(res$log2FoldChange),]
res_list <- res_down[rownames(toc_cutoff), ]$padj
names(res_list) <- rownames(res_down[rownames(toc_cutoff), ])
res_list <- res_list[!is.na(res_list)]
res_list_BP <- GOenrichment("KO_vs_WT_downregulated", geneList = res_list, ontology = "BP", cutoff = 0.05, nodeSize = 10, annType = annFUN.org, annObj = "org.Hs.eg.db", path = PATH_results, sig_nodes=5)

pdf(file = paste(PATH_results, "/KO_vs_WT_downregulated", "_", "BP", "_sigGO_graph",".pdf", sep=""))
showSigOfNodes(GOdata=res_list_BP$godata, termsP.value = score(res_list_BP$weightFisher), plotFunction = GOplot_custom, firstSigNodes = 5, useInfo = 'all', .NO.CHAR = 47, sigForAll=TRUE, sigFontSize = 70, sigBoxSize = c(36,31,31,14,14), elseFontSize = 18)
dev.off()

########################################################
##### DE genes annotated with interesting GO terms #####
########################################################

sensory_perception_pain <- c("GO:0051930")
potassium_ion_transport <- c("GO:0071805")
axonogenesis <- c("GO:0007409")
response_to_pain <- c("GO:0048265")
positive_regulation_synapse_assembly <- c("GO:0051965")

library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

sensory_perception_pain_up <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(sensory_perception_pain, rownames(res_up[res_up$padj < 0.05 & !is.na(res_up$padj),])), mart = ensembl)

potassium_ion_transport_up <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(potassium_ion_transport, rownames(res_up[res_up$padj < 0.05 & !is.na(res_up$padj),])), mart = ensembl)

axonogenesis_up <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(axonogenesis, rownames(res_up[res_up$padj < 0.05 & !is.na(res_up$padj),])), mart = ensembl)

response_to_pain_up <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(response_to_pain, rownames(res_up[res_up$padj < 0.05 & !is.na(res_up$padj),])), mart = ensembl)

positive_regulation_synapse_assembly_up <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(positive_regulation_synapse_assembly, rownames(res_up[res_up$padj < 0.05 & !is.na(res_up$padj),])), mart = ensembl)


sensory_perception_pain_down <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(sensory_perception_pain, rownames(res_down[res_down$padj < 0.05 & !is.na(res_down$padj),])), mart = ensembl)

potassium_ion_transport_down <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(potassium_ion_transport, rownames(res_down[res_down$padj < 0.05 & !is.na(res_down$padj),])), mart = ensembl)

axonogenesis_down <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(axonogenesis, rownames(res_down[res_down$padj < 0.05 & !is.na(res_down$padj),])), mart = ensembl)

response_to_pain_down <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(response_to_pain, rownames(res_down[res_down$padj < 0.05 & !is.na(res_down$padj),])), mart = ensembl)

positive_regulation_synapse_assembly_down <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), filters = c('go', 'ensembl_gene_id'), values = list(positive_regulation_synapse_assembly, rownames(res_down[res_down$padj < 0.05 & !is.na(res_down$padj),])), mart = ensembl)

write.csv(file = "sensory_perception_pain_up.csv", sensory_perception_pain_up)
write.csv(file = "potassium_ion_transport_up.csv", potassium_ion_transport_up)
write.csv(file = "axonogenesis_up.csv", axonogenesis_up)
write.csv(file = "response_to_pain_up.csv", response_to_pain_up)
write.csv(file = "positive_regulation_synapse_assembly_up.csv", positive_regulation_synapse_assembly_up)

write.csv(file = "sensory_perception_pain_down.csv", sensory_perception_pain_down)
write.csv(file = "potassium_ion_transport_down.csv", potassium_ion_transport_down)
write.csv(file = "axonogenesis_down.csv", axonogenesis_down)
write.csv(file = "response_to_pain_down.csv", response_to_pain_up)
write.csv(file = "positive_regulation_synapse_assembly_down.csv", positive_regulation_synapse_assembly_up)

###########
## WGCNA ##
###########
system("mkdir /home/george/Desktop/arc_ndcn_data/RNAseq/Network") 

ffun=filterfun(pOverA(p = 0.25, A = 10))
filt=genefilter(toc,ffun)
toc_coff=toc[filt,]

dds_coff <- DESeqDataSetFromMatrix(countData = toc_coff, colData = colData, design = ~ condition)

dds_coff$condition <- relevel(dds_coff$condition,"WT")

dds_coff <- estimateSizeFactors(dds_coff)
dds_coff <- estimateDispersions(dds_coff)

vsd_coff <- varianceStabilizingTransformation(dds_coff, blind=FALSE)
vsd_mat_coff <- assay(vsd_coff)
colnames(vsd_mat_coff) <- names(toc_coff)

save(vsd_mat_coff, vsd_coff, file = paste(PATH,"/Network/norm_counts_coff.RData", sep=""))

data_matrix <- vsd_mat_coff

WGCNA_matrix <- t(data_matrix[order(apply(data_matrix,1,mad), decreasing=T)[1:(nrow(data_matrix)/4)],])

s = abs(bicor(WGCNA_matrix))

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)

pdf(file = paste(PATH,"/Network/pick_soft_threshold_b.pdf", sep=""))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
dev.off()

gc()

#################
# Beta = lowest number for R2 0.9, find form diagram
beta <- 5
a <- s^beta
w <- 1-a

#create gene tree by average linkage hierarchical clustering 
geneTree = hclust(as.dist(w), method = 'average')

#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = w, deepSplit = 2, pamRespectsDendro = FALSE, cutHeight = 0.995, minClusterSize = 30)

#assign module colours
module.colours = labels2colors(modules)


library(ape)
#calculate eigengenes
MEList = moduleEigengenes(WGCNA_matrix, colors = module.colours)
MEs = MEList$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average')

MEDissThres = 0.2
# Plot the result
pdf(file = paste(PATH,"/Network/unmerged_modules_clustering.pdf", sep=""))
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

save(MEs, MEList, METree, modules, module.colours, geneTree, file = paste(PATH, "/Network/WGCNA.RData", sep=""))

#Merging modules
merge = mergeCloseModules(WGCNA_matrix, module.colours, cutHeight = MEDissThres, verbose = 3)
merged.colours = merge$colors;
mergedMEs = merge$newMEs;

pdf(file = paste(PATH, "/Network/Gene_dendro_merged_modules.pdf", sep=""))
plotDendroAndColors(geneTree, cbind(module.colours, merged.colours),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

module.colours = merged.colours
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(100));
moduleLabels = match(module.colours, colorOrder)-1;
MEs = mergedMEs;

module_annotation <- data.frame(moduleLabels = unique(moduleLabels), module.colours = unique(module.colours))

save(MEs, moduleLabels, module.colours, geneTree, module_annotation, file = paste(PATH, "/Network/merged_modules.RData", sep=""))

#reference genes = all ENSEMBL top 1/4 MAD genes 
ref_genes <- colnames(WGCNA_matrix)
ref_matrix <- WGCNA_matrix
ref_labels <- moduleLabels
ref_colours <- module.colours

#create data frame for GO analysis
library(org.Hs.eg.db)
GO = toTable(org.Hs.egGO); ENSEMBL = toTable(org.Hs.egENSEMBL)
GO_data_frame = data.frame(GO$go_id, GO$Evidence, ENSEMBL$ensembl_id[match(GO$gene_id,ENSEMBL$gene_id)])

GO_data_frame <- na.omit(GO_data_frame)

#create GOAllFrame object
library(AnnotationDbi)
GO_ALLFrame = GOAllFrame(GOFrame(GO_data_frame, organism = 'Homo Sapiens'))

#create gene set
library(GSEABase)
gsc <- GeneSetCollection(GO_ALLFrame, setType = GOCollection())

#perform GO enrichment analysis and save results to list - this make take several minutes
library(GOstats)
GSEAGO = vector('list',length(unique(ref_labels)))
for(i in 1:length(unique(ref_labels))){
  GSEAGO[[i]] = summary(hyperGTest(GSEAGOHyperGParams(name = 'Homo Sapiens GO', 
              geneSetCollection = gsc, geneIds = colnames(ref_matrix)[ref_labels==unique(ref_labels)[[i]]], 
              universeGeneIds = ref_genes, ontology = 'BP', pvalueCutoff = 0.05, 
              conditional = FALSE, testDirection = 'over')))
  print(i)
print(paste("Module:", unique(ref_labels)[[i]]))
}

ref_module_annotation <- data.frame(moduleLabels = unique(ref_labels), module.colours = unique(ref_colours))

cutoff_size = 100

GO_module_name = rep(NA,length(unique(unique(ref_labels))))
for (i in 1:length(unique(unique(ref_labels)))){
  GO_module_name[i] = 
    GSEAGO[[i]][GSEAGO[[i]]$Size<cutoff_size,
    ][which(GSEAGO[[i]][GSEAGO[[i]]$Size<cutoff_size,]$Count==max(GSEAGO[[i]][GSEAGO[[i]]$
    Size<cutoff_size,]$Count)),7]
}

module_annotation <- ref_module_annotation

module_annotation$GO_module  <- GO_module_name

save(gsc, GO_module_name, GSEAGO, module_annotation, file = paste(PATH,"/Network/GO_modules.RData", sep=""))

module_annotation <- module_annotation[order(module_annotation$moduleLabels), ]

rownames(module_annotation) <- module_annotation$moduleLabels

write.csv(module_annotation, file = paste(PATH,"/Network/module_annotation.csv", sep=""))

Genes_GO_module <- data.frame(LncRNA = colnames(WGCNA_matrix), GO_module = module_annotation[as.character(ref_labels),]$GO_module, module.label = ref_labels, module.colour = ref_colours)

Genes_GO_module$module.membership_weight <- NA

MM <- abs(bicor(WGCNA_matrix, MEs))

save(MM, module_annotation, file = paste(PATH,"/Network/Genes_ModuleMembership.RData", sep=""))

for(i in 1:nrow(Genes_GO_module)) {

Genes_GO_module$module.membership_weight[i] <- MM[Genes_GO_module$LncRNA[i], paste0("ME",Genes_GO_module$module.colour[i])]

}

write.csv(Genes_GO_module, file = paste(PATH,"/Network/Genes_GO_module.csv", sep=""))





