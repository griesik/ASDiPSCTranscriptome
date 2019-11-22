library('variancePartition')
library('edgeR')
library('BiocParallel')
library("mgsa")
library("biomaRt")
library(dplyr)

setwd("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA")
Neuron=read.table("countData_RUV_Normalized_Neurons_nDEGs_HKgenes_p0.4_k4.txt")
sampleData=read.delim("sample_sheet_all_cells.txt")
sampleData=sampleData[c(40, 41, 44, 45, 46, 48:54, 56:62),]
rownames(sampleData)=sampleData[,1]
sampleData=sampleData[,-1]
geneExpr=DGEList(Neuron)
geneExpr = calcNormFactors( geneExpr )

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)
# The variable to be tested must be a fixed effect
form= ~ condition + (1|sample)
# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, sampleData )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, sampleData )
# Examine design matrix
head(fitmm$design)
# Get results of hypothesis test on coefficients of interest

topTable( fitmm, coef='conditionP')
Neuronresult=as.data.frame(topTable (fitmm, coef="conditionP", number=15026))

#annotating the table
ensembl <- rownames(Neuronresult)
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
geneR=geneR[,c(2,1,3)]
Neuronresult[,7]=rownames(Neuronresult)
colnames(Neuronresult)[7]="Ensembl"
Neuronresult <- left_join (Neuronresult, geneR, by = "Ensembl")
Neuronresult=Neuronresult[,c(7,8,1:6)]
write.table(Neuronresult, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/Neurons/Neuron_dream_result.txt")



#####limma analysis####
# apply duplicateCorrelation is two rounds
design = model.matrix( ~ condition, sampleData)
vobj_tmp = voom( geneExpr, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=sampleData$sample)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
vobj = voom( geneExpr, design, plot=FALSE, block=sampleData$sample, correlation=dupcor$consensus)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor <- duplicateCorrelation(vobj, design, block=sampleData$sample)

# But this step uses only the genome-wide average for the random effect
fitDupCor <- lmFit(vobj, design, block=sampleData$sample, correlation=dupcor$consensus)

# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes( fitDupCor )
Neuronresultdupcor=toptable(fitDupCor, coef="conditionP", number=15026)

#annotating the table
ensembl <- rownames(Neuronresultdupcor)
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
geneR=geneR[,c(2,1,3)]
Neuronresultdupcor[,6]=rownames(Neuronresultdupcor)
colnames(Neuronresultdupcor)[6]="Ensembl"
Neuronresultdupcor <- left_join (Neuronresultdupcor, geneR, by = "Ensembl")
Neuronresultdupcor=Neuronresultdupcor[,c(6,7,1:5)]
write.table(Neuronresultdupcor, file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/Neurons/Neuron_dupcor_result.txt")




library(mgsa)
goterms=readGAF("goa_human.gaf", evidence=NULL, aspect="P")

#preparing the dream result table for enrichment analysis
#Neuronresult=read.delim("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/Neurons/Neuron_dream_result.txt", sep="")
Neuronresult[,7]=rownames(Neuronresult)
colnames(Neuronresult)[7]="Ensembl"
Neuronresult=Neuronresult[,c(7,1:6)]
ensembl <- Neuronresult[,1]
geneR <-if (interactive()) {
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"),
        filters    = "ensembl_gene_id",
        values     = ensembl , 
        mart       = mart)
}

colnames (geneR) <- c("Gene", "Ensembl", "Entrez")
geneR=geneR[,c(2,1,3)]

logtrans <- inner_join (Neuronresult, geneR, by = "Ensembl")
logtrans=logtrans[!is.na(logtrans$Entrez),]
logtrans=distinct(logtrans,Entrez, .keep_all= TRUE)
logtrans=logtrans[,c(9,1,8,2:7)]
up=which(logtrans[,6]>=0)
dreamup=logtrans[up,]
#take the colum with t statistics
statisticup=dreamup[,6]
names(statisticup)=dreamup[,1]
down=which(logtrans[,6]<=0)
dreamdown=logtrans[down,]
#take the colum with t statistics
statisticdown=dreamdown[,6]
names(statisticdown)=dreamdown[,1]




sets=goterms@sets

for(j in c(1:16114))
{sets[[j]]= subset (sets[[j]], sets[[j]] %in% logtrans[,1])
  }


enrichdown=cameraPR(statisticdown, sets, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)
enrichup=cameraPR(statisticup, sets, use.ranks = TRUE, inter.gene.cor=0.01, sort = FALSE)




###preparing the dupCorr results table for enrichment analysis
ensembl <- rownames(Neuronresultdupcor)
geneR <-if (interactive()) {
      mart <- useMart("ensembl")
       datasets <- listDatasets(mart)
       mart <- useDataset("hsapiens_gene_ensembl",mart)
     getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"),
                         filters    = "ensembl_gene_id",
                         values     = ensembl , 
                       mart       = mart)
   }

colnames (geneR) <- c("Gene", "Ensembl", "Entrez")
geneR=geneR[,c(2,1,3)]
Neuronresultdupcor[,6]=rownames(Neuronresultdupcor)
colnames(Neuronresultdupcor)[6]="Ensembl"
logtrans <- left_join (Neuronresultdupcor, geneR, by = "Ensembl")
logtrans=distinct(logtrans,Entrez, .keep_all= TRUE)
statistic=logtrans[,2]
names(statistic)=logtrans[,8]



sets=goterms@sets
enrichdown=cameraPR(statisticdown, sets, use.ranks = TRUE, inter.gene.cor=0.01, sort = TRUE)


#####rankprod anakysis####
#Neurons

Neuronlog=read.csv("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/Neurons/logtransfcounts_notCol_Neurons_nDEGs_HKgenes_p0.4_k4.csv")
sampleData=read.delim("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/DATA/sample_sheet_all_cells.txt")
sampleData <- sampleData [c(40, 41, 44, 45, 46, 48:54, 56:62),]
genes=Neuronlog[,2]
Neuronlog=Neuronlog[,c(10:28)]
class=sampleData[,5]
RP.out <- RankProducts(Neuronlog,class, logged=TRUE, na.rm=FALSE,plot=FALSE, rand=123)
result=topGene(RP.out,cutoff=0.05,method="pfp",logged=TRUE,logbase=2,gene.names=genes)
resultdown=result$Table1
resultup=result$Table2
