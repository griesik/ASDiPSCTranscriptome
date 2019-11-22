####NPC Analysis - FINAL SCRIPT ###########
## Load required packages ##########################################
library(edgeR)
library(RUVSeq)
library(EDASeq)
library(WGCNA)
library(calibrate)
library(DESeq2)
library(pheatmap)
library (RColorBrewer)
library (dplyr)
library (biomaRt)
library(multtest)
library(gplots)
library(janitor)
library(variancePartition)
library(BiocParallel)
library(fgsea)



####Loading Data######
#Load FPKM table to retain only expressed genes (FPKM>1 in more than half of the samples, calculated seperately to NPC and neurons)
setwd ("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA")
geneRPKM=read.table ("allData_FPKM_renormalized_IV.txt", sep= "\t", header=TRUE)
rownames(geneRPKM)=geneRPKM [,1]
geneRPKM = geneRPKM [,-1]

#Load list of HK genes (house keeping genes selected from Eisenberg & Levanon, 2013)
HKgenes_full <- read.csv ("HK_full_gene_list_ensembl_biomart.csv")

#load sampleData
sampleInfo <- read.table ("sample_sheet_all_cells.txt", header = TRUE)

#load NPC countdata
counts=read.table("countdata_20M_NPC.txt", header=TRUE)


#exclude samples with more than 55% in neuronal proportion

countData <- counts [,c(1:10,12:14,22:37)]
sampleData <- sampleInfo [c(1:10,12:14,22:37),]

#filtering genes with FPKM<1 in more than half of the samples
geneRPKM2 <- geneRPKM [,c(1:10,12:14,22:37)]
filtered <- which(rowSums(geneRPKM2 > 1) >= 15)
RPKM <- geneRPKM2[filtered,]
countData=subset (countData, rownames(countData) %in% rownames(RPKM))



####Plot PCA - before normalization#####
#SuppFig1
m=match(colnames(countData), sampleData$label)
pdf("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/BeforeNorm_NPC.pdf")
boxplot(log2(countData))
pca=princomp(cor(countData, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,BeforeNorm, batch", pch=20, col=labels2colors(sampleData$batch[m]))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,BeforeNorm, Neur_proportion", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("batch"))
plotPCA(vsd, intgroup=c("Neur_prop_percentage"))
plotPCA(vsd, intgroup=c("condition"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
sampleTrees = hclust(dist(t(as.data.frame(assay(vsd)))), method = "average")
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
plot(sampleTrees, main = paste("Sample clustering NPC"),
     xlab="", sub="", cex = 0.7)
dev.off()


####RUVseq normalization#####
#find the intersect between countData and HKgenes
HKgenes <- intersect(HKgenes_full$Ensembl_biomart, rownames(countData))

#run a first differential expression analysis to find non-DEGs in the non-normalized dataset
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)
dds$condition <- relevel(dds$condition, ref="C")
ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq)
summary(res)
res <- as.data.frame(res)

####set nonDEGs#####
#change the pvalue from 0.1 to 0.9, by 0.1, to test the best set of genes for normalization
nDEGs_deseq <- subset (res, res$pvalue>=0.3)
nDEGs_deseq <- rownames(nDEGs_deseq)

#intersect nonDEGs with HKgenes (trying to keep only those house keep genes that have a consistent expression)
controlGenes <- intersect(nDEGs_deseq,HKgenes)


m=match(colnames(countData), sampleData$label)

#####Run RUVseq analysis and plot results####
pdf("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/nDEGs_HKgenes_p0.3_NPC.pdf")

RUV1 <- RUVg(as.matrix(countData), controlGenes, k=1)
N1 <- RUV1$normalizedCounts

boxplot(log2(N1))
pca=princomp(cor(N1, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k1, batch", pch=20, col=labels2colors(sampleData$batch[m]))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k1, Neur_prop", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N1, colData = sampleData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("batch"))
plotPCA(vsd, intgroup=c("Neur_prop_percentage"))
plotPCA(vsd, intgroup=c("condition"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
sampleTrees = hclust(dist(t(as.data.frame(assay(vsd)))), method = "average")
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
plot(sampleTrees, main = paste("Sample clustering NPC"),
     xlab="", sub="", cex = 0.7)



#Supp Fig.3
RUV2 <- RUVg(as.matrix(countData), controlGenes, k=2)
N2 <- RUV2$normalizedCounts

boxplot(log2(N2))
pca=princomp(cor(N2, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k2, batch", pch=20, col=labels2colors(sampleData$batch[m]))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k2, Neur_prop_NPC", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N2, colData = sampleData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("batch"))
plotPCA(vsd, intgroup=c("Neur_prop_percentage"))
plotPCA(vsd, intgroup=c("condition"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
sampleTrees = hclust(dist(t(as.data.frame(assay(vsd)))), method = "average")
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
plot(sampleTrees, main = paste("Sample clustering NPC"),
     xlab="", sub="", cex = 0.7)


RUV3 <- RUVg(as.matrix(countData), controlGenes, k=3)
N3 <- RUV3$normalizedCounts

boxplot(log2(N3))
pca=princomp(cor(N3, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k3, batch", pch=20, col=labels2colors(sampleData$batch[m]))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k3, Neur_prop", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N3, colData = sampleData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("batch"))
plotPCA(vsd, intgroup=c("Neur_prop_percentage"))
plotPCA(vsd, intgroup=c("condition"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
sampleTrees = hclust(dist(t(as.data.frame(assay(vsd)))), method = "average")
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
plot(sampleTrees, main = paste("Sample clustering NPC"),
     xlab="", sub="", cex = 0.7)



RUV4 <- RUVg(as.matrix(countData), controlGenes, k=4)
N4 <- RUV4$normalizedCounts

boxplot(log2(N4))
pca=princomp(cor(N4, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k4, batch", pch=20, col=labels2colors(sampleData$batch[m]))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k4, Neur_prop_NPC", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N4, colData = sampleData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("batch"))
plotPCA(vsd, intgroup=c("Neur_prop_percentage"))
plotPCA(vsd, intgroup=c("condition"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
sampleTrees = hclust(dist(t(as.data.frame(assay(vsd)))), method = "average")
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
plot(sampleTrees, main = paste("Sample clustering NPC"),
     xlab="", sub="", cex = 0.7)



RUV5 <- RUVg(as.matrix(countData), controlGenes, k=5)
N5 <- RUV5$normalizedCounts


boxplot(log2(N5))
pca=princomp(cor(N5, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k5, batch", pch=20, col=labels2colors(sampleData$batch[m]))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k5, cell", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N5, colData = sampleData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("batch"))
plotPCA(vsd, intgroup=c("Neur_prop_percentage"))
plotPCA(vsd, intgroup=c("condition"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
sampleTrees = hclust(dist(t(as.data.frame(assay(vsd)))), method = "average")
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
plot(sampleTrees, main = paste("Sample clustering NPC"),
     xlab="", sub="", cex = 0.7)
dev.off()

####calculating correlation to neuronal proportion####
#SuppTable2
pcavalues=c(1:5)
pcavalues=as.data.frame(pcavalues)
pcavalues$pca1cor = "NA"
pcavalues$pca1pvalue = "NA"
pcavalues$pca2cor = "NA"
pcavalues$pca2pvalue = "NA"
NPC_prop=sampleData[,10]

#set pvalue, varying from 0.1 to 0.9, by=0.1
nDEGs_deseq <- subset (res, res$pvalue>=0.1)
nDEGs_deseq <- rownames(nDEGs_deseq)
controlGenes <- intersect(nDEGs_deseq,HKgenes)


for (j in c(1:nrow(pcavalues)))
{RUV <- RUVg(as.matrix(countData), controlGenes, k=j)
N <- RUV$normalizedCounts
pca=princomp(cor(N, method="s"))
cor1=cor(pca$loadings[,1], NPC_prop)
cor2=cor(pca$loadings[,2], NPC_prop)
p1=corPvalueStudent(cor1,29)
p2=corPvalueStudent(cor2,29)
pcavalues$pca1cor[j] = cor1
pcavalues$pca1pvalue[j] = p1
pcavalues$pca2cor[j] = cor2
pcavalues$pca2pvalue[j] = p2}

#input these results to a new object for every pvalue tested
pcavalues1=pcavalues

pcavaluesfinal=cbind(pcavalues,pcavalues1,pcavalues2,pcavalues3,pcavalues4,pcavalues5,pcavalues6,pcavalues7,pcavalues8,pcavalues9)
write.csv(pcavaluesfinal, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/pca_correlation_values_nDEGs_HKgenes_NPC.csv")

#calculate the pca values for the non-normalized data
pca=princomp(cor(countData, method="s"))
cor1=cor(pca$loadings[,1], Neur_prop)
cor2=cor(pca$loadings[,2], Neur_prop)
p1=corPvalueStudent(cor1,29)
p2=corPvalueStudent(cor2,29)
print(cor1)
print(cor2)
print(p1)
print(p2)



####dream analysis####
#after chosing the best parameters based on RUVseq analysis, use the normalized count Data
#generated by these parameters
#In case of NPCs, the best parameters were: pvalue to choose the non-DEGs:0.3; k=2

nDEGs_deseq <- subset (res, res$pvalue>=0.3)
nDEGs_deseq <- rownames(nDEGs_deseq)

controlGenes <- intersect(nDEGs_deseq,HKgenes)

RUV <- RUVg(as.matrix(countData), controlGenes, k=2)
N <- RUV$normalizedCounts

write.table (N, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/countData_RUV_Normalized_NPC_nDEGs_HKgenes_p0.3_k2.txt")

geneExpr=DGEList(N)
geneExpr = calcNormFactors( geneExpr )
rownames(sampleData)=sampleData[,1]
sampleData1=sampleData[,-1]

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)
# The variable to be tested must be a fixed effect
form= ~ condition + (1|sample)
# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, sampleData1 )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, sampleData1 )
# Examine design matrix
head(fitmm$design)
# Get results of hypothesis test on coefficients of interest

NPCresult=as.data.frame(topTable (fitmm, coef="conditionP", number=13818))


NPCresult[,7]=rownames(NPCresult)
colnames(NPCresult)[7]="Ensembl"
NPCresult=NPCresult[,c(7,1:6)]
ensembl <- NPCresult[,1]
geneR <-if (interactive()) {
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype"),
        filters    = "ensembl_gene_id",
        values     = ensembl , 
        mart       = mart)
}

colnames (geneR) <- c("Gene", "Ensembl", "Gene_type")



#Supp Table 6
NPCresult <- left_join (NPCresult, geneR, by = "Ensembl")
NPCresult[,c(10,11)]="NA"
colnames(NPCresult)[c(10,11)]=c("rawFC","FoldChange")
for(j in c(1:nrow(NPCresult)))
{NPCresult$rawFC[j]= 2^NPCresult$logFC[j]}
NPCresult[,10]=as.numeric(NPCresult$rawFC)
for(j in c(1:nrow(NPCresult)))
{NPCresult$FoldChange[j]= if (NPCresult$rawFC[j]>1) {NPCresult$rawFC[j]} else {-1/NPCresult$rawFC[j]}}

#load Sfari and ID genes data
sfari=read.csv("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/SFARI-Gene_genes_11-01-2019_ens.csv")
sfari=sfari[,c(1,6)]
colnames(sfari)[2]="Sfari_score"
ID=read.csv("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/ID_genes_2019.csv")
ID=ID[,c(1,6)]
colnames(ID)[2]="ID_score"
NPCresult=left_join(NPCresult, sfari)
NPCresult=left_join(NPCresult,ID)
NPCresult=NPCresult[,c(1,8,9,2,10,11,3:7,12,13)]

write.csv( as.data.frame(NPCresult), file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/dreamResults_NPC_nDEGs_HKgenes_p0.3_k2.csv", row.names = FALSE )



#####enrichment analysis####

#preparing the dream result table for enrichment analysis

#creating columns with entrez ID and gene type

geneR <-if (interactive()) {
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
        filters    = "ensembl_gene_id",
        values     = ensembl , 
        mart       = mart)
}
colnames (geneR) <- c("Ensembl", "Entrez")
enrTable <- inner_join (NPCresult, geneR, by = "Ensembl")
enrTable=enrTable[!is.na(enrTable$Entrez),]
enrTable=distinct(enrTable,Entrez, .keep_all= TRUE)
enrTable=enrTable[,c(13,1:12)]
#subset up and down regulated genes 
up=which(enrTable[,8]>=0)
dreamup=enrTable[up,]
#take the colum with t statistics
statisticup=dreamup[,8]
names(statisticup)=dreamup[,1]
down=which(enrTable[,8]<=0)
dreamdown=enrTable[down,]
#take the colum with t statistics
statisticdown=dreamdown[,8]
names(statisticdown)=dreamdown[,1]

#msigDB lists
#pathways
index=gmtPathways("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/msigDB_pathways.gmt")

for(j in c(1:2199))
{index[[j]]= subset (index[[j]], index[[j]] %in% enrTable[,1])
}

msigPathwaysDown=cameraPR(statisticdown, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)
msigPathwaysUp=cameraPR(statisticup, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)

write.csv(msigPathwaysDown, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_down_pathways.csv")
write.csv(msigPathwaysUp, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_up_pathways.csv")


#GO_BP
index=gmtPathways("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/msigDB_GO_BP.gmt")

for(j in c(1:7350))
{index[[j]]= subset (index[[j]], index[[j]] %in% enrTable[,1])
}

msigGOBPDown=cameraPR(statisticdown, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)
msigGOBPUp=cameraPR(statisticup, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)

write.csv(msigGOBPDown, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_down_GOBP.csv")
write.csv(msigGOBPUp, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_up_GOBP.csv")


#microRNA_targets
index=gmtPathways("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/msigDB_microRNA_targets.gmt")

for(j in c(1:221))
{index[[j]]= subset (index[[j]], index[[j]] %in% enrTable[,1])
}

msigmicroRNADown=cameraPR(statisticdown, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)
msigmicroRNAUp=cameraPR(statisticup, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)

write.csv(msigmicroRNADown, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_down_microRNA.csv")
write.csv(msigmicroRNAUp, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_up_microRNA.csv")


#TF motifs
index=gmtPathways("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/msigDB_TF_motifs.gmt")

for(j in c(1:610))
{index[[j]]= subset (index[[j]], index[[j]] %in% enrTable[,1])
}

msigTFDown=cameraPR(statisticdown, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)
msigTFUp=cameraPR(statisticup, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)

write.csv(msigTFDown, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_down_TFmotifs.csv")
write.csv(msigTFUp, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_up_TFmotifs.csv")


#hallmark
index=gmtPathways("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/msigDB_hallmark.gmt")

for(j in c(1:50))
{index[[j]]= subset (index[[j]], index[[j]] %in% enrTable[,1])
}

msigHallmarkDown=cameraPR(statisticdown, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)
msigHallmarkUp=cameraPR(statisticup, index, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)

write.csv(msigHallmarkDown, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_down_hallmark.csv")
write.csv(msigHallmarkUp, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/enrichment_NPC_up_hallmark.csv")


#####generate log transformed counts####
#create deseq object
dds <- DESeqDataSetFromMatrix(countData = N, colData = sampleData, design = ~ condition)

dds$condition <- relevel(dds$condition, ref="C")

#to extract the log-transformed counts
vsd<- varianceStabilizingTransformation(dds)
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)
ncvsd2<- ncvsd
ncvsd2[,30]<-rownames (ncvsd2)
colnames(ncvsd2)[30]="Ensembl"
logtrans <- left_join(NPCresult, ncvsd2, by = "Ensembl")
write.csv(logtrans, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/logtransfcounts_NPC_nDEGs_HKgenes_p0.3_k2.csv", row.names=FALSE)


#####WGCNA analysis######
#check for variance over the mean and possibliy remove some genes
variance <- apply(ncvsd,1,sd)/apply(ncvsd,1,mean)
hist(variance)

keep=which(apply(ncvsd,1,sd)/apply(ncvsd,1,mean)>= 0.0175)#I calculated the variance over the mean and tried different cuttoffs to see how many genes were still retained
ncvsd=ncvsd[keep,]


#transpose the table
dat=t(ncvsd)
infoData=rownames(ncvsd)



#Run WGCNA#
# run soft pick threshold
# Test a series of powers to which co-expression similarity is raised to calculate adjacency
powers=c(c(1:10), seq(from = 12, to = 20, by = 2)) 
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5, blockSize = 20000, networkType = "signed") # Signed gives an indication of positive vs negative correlations
#Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9 #what is that?
# Fit to scale-free topology 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence")); 
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = cex1, col="red"); 
abline(h=c(0.8, 0.5), col="red")
# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold (power)", 
     ylab = "Mean Connectivity", type = "n", main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col="red")
pdf(file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/WGCNA_NPC_SoftThreshold_nDEGs_HKgenes_p0.3_k2.pdf")
dev.off()

#construction of the network
# Chosen soft threshold in this case is 12; it sets the fit to Scale-free Topology > 0.8 while retaining as much connectivity as possible
net=blockwiseModules(dat, power=10, numericLabels=TRUE, networkType = "signed", 
                     minModuleSize=50, mergeCutHeight=0.15, saveTOMs=FALSE, verbose=6,minKMEtoStay = 0.5,
                     nThreads=24, maxBlockSize=20000, checkMissingData=FALSE)


## Labelling the modules with a colour tag
modules=as.data.frame(table(net$colors)); colnames(modules)=c("Label", "N") 
modules$Label=paste("M", modules$Label, sep=""); 
modules$Color=c("grey",labels2colors(modules$Label[-1])) 
moduleLabel=paste("M",net$colors, sep=""); 
moduleColor=modules$Color[match(moduleLabel, modules$Label)] 


####dendogram modules
sizeGrWindow(12, 9)
pdf(file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/Cluster_dendogram_NPC.pdf")

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#####get kMe table####
## Calculating kMEs
KMEs<-signedKME(dat, net$MEs,outputColumnName = "M") # signedKME calculates eigengene-based connectivity i.e. module membership. Also, outputColumnName gives us a prefix. Also, the corFnc defaults to Pearson
kme=data.frame(infoData[match(colnames(dat), infoData)], moduleColor,moduleLabel, KMEs)
colnames(kme)[1]="Ensembl"


#Re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue <0.05, and kME>0.5, will be moved to a junk module
kmeInfoCols=c(1:3)
kmedata=kme[,-kmeInfoCols]; 
pvalBH=kmedata; pvalBH[,]=NA 
for (j in c(1:ncol(pvalBH)))
{
  p=mt.rawp2adjp(corPvalueStudent(kmedata[,j], nSamples=29), proc="BH")
  pvalBH[,j]=p$adjp[order(p$index),2]
}

kme$newModule="NA" 
for (j in c(1:nrow(kmedata)) )
{
  if (j==1) print("Working on genes 1:10000"); if(j==10000) print(paste("Working on genes 10000:", nrow(kmedata)))
  m=which(kmedata[j,]==max(kmedata[j,]))
  if ((pvalBH[j,m]<0.05)&(kmedata[j,m]>0.5)) kme$newModule[j]=as.character(colnames(kmedata)[m])
}

## Assign genes not associated to any module to M0
kme$newModule[which(kme$newModule%in%"NA")]="M0" 

## Replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
kme$newColor=kme$moduleColor[match(kme$newModule, kme$moduleLabel)] 
kme$moduleLabel=kme$newModule; kme$moduleColor=kme$newColor 
kme=kme[,-grep("newModule", colnames(kme))];kme=kme[,-grep("newColor", colnames(kme))] 

#Saving kMEs 
mod=modules$Label[-1] 
kmeTable=kme[,kmeInfoCols]; 
for(j in c(1:length(mod)))
{
  kmeTable=cbind(kmeTable, kmedata[,match(mod[j],colnames(kmedata))]);colnames(kmeTable)[ncol(kmeTable)]=paste("kME", mod[j], sep="_")                                                                             
  kmeTable=cbind(kmeTable, pvalBH[,match(mod[j],colnames(pvalBH))]);colnames(kmeTable)[ncol(kmeTable)]=paste("pvalBH", mod[j], sep="_")
}

kmeTable=left_join(kmeTable, geneR, by="Ensembl")
kmeTable=kmeTable[,c(1,26,2:25)]
#Sup Table 9
write.csv(kmeTable, "C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/kME_NPC_nDEGs_HKgene_p0.3_k2.csv", row.names=FALSE)

#Saving Module Eigengenes 
me<-data.frame(rownames(dat), net$MEs) # Sample names bound to module eigengenes
colnames(me)[-1]=gsub("ME", "M", colnames(me)[-1])
colnames(me)[1]="Sample"
write.csv(me, "C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/ME_NPC_nDEGs_HKgene_p0.3_k2.csv", row.names=FALSE)

#Plotting Module Eigengene Values

colors=rep("turquoise", nrow(me) )
colors[grep("_P_", me$Sample)]="red"

pdf("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/moduleBarplots_NPC_nDEGs_HKgene_p0.3_k2.pdf", height=5, width=15)

mod=paste("M", c(0:(ncol(me)-2)), sep="")
for(m in mod) 
{ 
  j=match(m, colnames(me))
  col=kme$moduleColor[match(m, kme$moduleLabel)]
  barplot(me[,j],  xlab="Samples", ylab="ME",col=colors, main=m, names=me[,1], cex.names=0.5, axisnames = FALSE)
  
}
dev.off()

# get biological variables correlations
# Define numbers of genes and samples
Samples = rownames(dat);
traitRows = match(Samples, sampleData$label);
SampleInfo = sampleData[traitRows, -1];
rownames(SampleInfo) = sampleData[traitRows, 1];
nGenes = ncol(dat);
nSamples = nrow(dat);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dat, moduleColor)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, SampleInfo, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
pdf("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/Module-TRait_NPC_nDEG_HKgenes_p0.3_k2.pdf")
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(SampleInfo),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#####Binding modules names to the dream NPCresult table#####
#joining kmeTable with genes excluded in variance filter
modules_griesi=kmeTable[,c(1,3,4)]
NPCresult=left_join(NPCresult,modules_griesi)
NPCresult[["moduleLabel"]][is.na(NPCresult[["moduleLabel"]])] = "M0"
NPCresult[["moduleColor"]][is.na(NPCresult[["moduleColor"]])] = "grey"

#supp table 6
write.csv(as.data.frame(NPCresult),file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC/dreamResults_NPC_nDEGs_HKgenes_p0.3_k2.csv", row.names = FALSE)


####overlapping of modules with external data####

geneR=NPCresult[,1]
labelR=NPCresult[,15]

#setting the list of modules from the other studies
setwd("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA")
gupta=read.csv("modules_Gupta.csv")
gandal=read.csv("modules_Gandal.csv")
mariani=read.csv("modules_mariani_ensembl.csv")
tylee=read.csv("modules_tylee_ensembl.csv")
voineagu=read.csv("modules_voineagu_ensembl.csv")
schafer=read.csv("modules_schafer.csv")
parikshak2013=read.csv("modules_Parikshak2013.csv")
parikshak2016=read.csv("modules_parikshak2016.csv")
derosa35=read.csv("modules_DeRosa35DIV.csv")
#roussos=read.csv("modules_roussos.csv")
#derosa135=read.csv("modules_DeRosa135DIV.csv")
#lin=read.csv("modules_Lin.csv")
#lombardo=read.csv("modules_lombardo_ensembl.csv")

allStudies=rbind(gupta,gandal,mariani,tylee,voineagu,parikshak2013,parikshak2016,
                 derosa35)

write.csv(allStudies, file="modules_allStudies.csv", row.names = FALSE)


#running overlap test
setwd("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/NPC")


userListEnrichment(
  geneR, labelR, 
  fnIn = ("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/modules_allStudies.csv"),
  nameOut = "module_overlap_npc_allStudies_final.csv",
  omitCategories = c("M0","grey"), outputGenes = TRUE) 


#calculating odds ratio
#to calculate all the fields for the fisher table
otherstudy2=left_join(schafer,NPCresult, by="Gene")
ss=subset(otherstudy2, otherstudy2$Module.x=="saddlebrown_Schafer" & otherstudy2$Module.y=="M10")
sn=subset(otherstudy2, otherstudy2$Module.x=="saddlebrown_Schafer" & !otherstudy2$Module.y=="M10")
ns=subset(otherstudy2, !otherstudy2$Module.x=="saddlebrown_Schafer" & otherstudy2$Module.y=="M10")
nn=subset(otherstudy2, !otherstudy2$Module.x=="saddlebrown_Schafer" & !otherstudy2$Module.y=="M10")

enrichmenttab=read.csv("module_overlap_npc_allStudies_final.csv")
tabfisher=read.delim("tablefisher_module_overlap_npc.txt")
enrichmenttab=subset(enrichmenttab, enrichmenttab$InputCategories=="M10")
enrichmenttab=left_join(enrichmenttab,tabfisher, by="UserDefinedCategories")
enrichmenttab=enrichmenttab[!is.na(enrichmenttab$yesother.yesM1),]


enrichmenttab$oddsratio = "NA"
for(j in c(1:nrow(enrichmenttab)))
{counts = (matrix(data = as.numeric(as.character(enrichmenttab[j,c(8,9,10,11)])), nrow=2))
res = fisher.test(counts)
enrichmenttab$oddsratio[j]=res$estimate}
enrichmenttab=enrichmenttab[,c(1:4,12,5:11)]
write.csv(enrichmenttab, file="module_overlap_npc_allStudies_M10_final.csv")

