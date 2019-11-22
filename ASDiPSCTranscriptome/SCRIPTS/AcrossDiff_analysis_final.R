###Differentiation Analysis - FINAL SCRIPT######
#Analysis of differential expression between NPC and Neurons
##Load libraries####
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



####Differential expression analysis with normalized data#####
####Differential expression between NPC and neurons#####
####Loading Data
#Fig.1C
setwd("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/artigo_expressao/ASDiPSCTranscriptome")
NPCData <- read.table ("DATA/countData_RUV_Normalized_NPC_nDEGs_HKgenes_p0.3_k2.txt", header=TRUE)
NeuronData <- read.table ("DATA/countData_RUV_Normalized_Neurons_nDEGs_HKgenes_p0.4_k4.txt", header = TRUE)

NPCData[,30]=rownames(NPCData)
colnames(NPCData)[30]="Gene"
NeuronData[,20]=rownames(NeuronData)
colnames(NeuronData)[20]="Gene"
countData = full_join (NPCData,NeuronData, by="Gene")
rownames(countData)=countData[,30]
countData=countData[,-30]
countData[is.na(countData)]=0

write.table (countData, file="DATA/countData_RUV_normalized_allCells.txt")


sampleInfo <- read.table ("DATA/sample_sheet_all_cells.txt", header = TRUE)
sampleData <- sampleInfo [c(1:29,30, 31, 34, 35, 36, 38:44, 46:52),]

#####Dream analysis####
geneExpr=DGEList(countData)
geneExpr = calcNormFactors( geneExpr )
rownames(sampleData)=sampleData[,1]
sampleData1=sampleData[,-1]

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)
# The variable to be tested must be a fixed effect
form= ~ cell_num + (1|sample)
# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, sampleData1 )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, sampleData1 )
# Examine design matrix
head(fitmm$design)
# Get results of hypothesis test on coefficients of interest

result=as.data.frame(topTable (fitmm, coef="cell_num", number=16368))


result[,7]=rownames(result)
colnames(result)[7]="Ensembl"
result=result[,c(7,1:6)]
ensembl <- result[,1]
geneR <-if (interactive()) {
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
        filters    = "ensembl_gene_id",
        values     = ensembl , 
        mart       = mart)
}

colnames (geneR) <- c("Gene", "Ensembl")



#Supp Table 5
result <- left_join (result, geneR, by = "Ensembl")
result[,c(9,10)]="NA"
colnames(result)[c(9,10)]=c("rawFC","FoldChange")
for(j in c(1:nrow(result)))
{result$rawFC[j]= 2^result$logFC[j]}
result[,9]=as.numeric(result$rawFC)
for(j in c(1:nrow(result)))
{result$FoldChange[j]= if (result$rawFC[j]>1) {result$rawFC[j]} else {-1/result$rawFC[j]}}


result=result[,c(1,8,2,9,10,3:7)]

write.csv( as.data.frame(result), file="RESULTS/Across_diff/dreamanalysis_NPCxNeurons.csv", row.names = FALSE )



#####generate log transformed counts####
#create deseq object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ cell)

dds$condition <- relevel(dds$cell, ref="NPC")

#to extract the log-transformed counts
vsd<- varianceStabilizingTransformation(dds)
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)
ncvsd2<- ncvsd
ncvsd2[,49]<-rownames (ncvsd2)
colnames(ncvsd2)[49]="Ensembl"
logtrans <- left_join(result, ncvsd2, by = "Ensembl")
write.csv(logtrans, file="RESULTS/Across_diff/logtransfcounts_Across_diff_nDEGs_HKgenes.csv", row.names=FALSE)






###preparing data for CoNTExT analysis####

##finding gene symbols and entrez gene numbers
#the object ensembl has already been defined above
geneR <-if (interactive()) {
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"),
        filters    = "ensembl_gene_id",
        values     = ensembl , 
        mart       = mart)
}
ncvsd[,49]=rownames(ncvsd)
colnames(ncvsd)[49]="ensembl_gene_id"


####
Contmatrix=inner_join(ncvsd, geneR, by="ensembl_gene_id")
#removing entries with no entrez ID
Contmatrix=Contmatrix[!is.na(Contmatrix$entrezgene_id),]
#remove repeated entrez ID
Contmatrix=distinct(Contmatrix,entrezgene_id, .keep_all= TRUE)


rownames(Contmatrix)=Contmatrix[,51]
annot=Contmatrix[,c(51,50,51)]
rownames(annot)=annot[,1]
annot=annot[,-1]
colnames(annot)=c("SYMBOL","ENTREZ_GENE_ID")
Contmatrix=Contmatrix[,-c(49,50,51)]

write.csv(Contmatrix, file="RESULTS/Across_diff/logData_allcells_onlyEntrez_context.txt", row.names=FALSE)
write.csv(annot,file="RESULTS/Across_diff/annotation_file_context_entrez.txt", row.names=FALSE)


#####cell marker characterization####
panglao=read.delim("DATA/PanglaoDB_cell_markers_filtered_list.txt")
#ensembl=panglao[,2]
#mart <- useMart("ensembl")
#datasets <- listDatasets(mart)
#mart <- useDataset("hsapiens_gene_ensembl",mart)
#panglaoconv=  getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#  filters    = "hgnc_symbol",
# values     = ensembl , 
#  mart       = mart)


logData=read.csv("RESULTS/Across_diff/logtransfcounts_Across_diff_nDEGs_HKgenes.csv")

panglao_expdata=left_join(panglao,logData, by="Ensembl")
write.table(panglao_expdata, file="RESULTS/Across_diff/cell_markers_exp_data.txt")
