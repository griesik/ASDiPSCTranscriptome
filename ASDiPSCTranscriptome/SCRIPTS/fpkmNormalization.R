####FPKM Calculation#####
#This script is to normalize count numbers by exonic lenght and library size, generating FPKM measures

#Load libraries
library (dplyr)
library(data.table)

#combining NPC and Neuron data
setwd("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/artigo_expressao/ASDiPSCTranscriptome")
countData = read.table("DATA/countdata_full.txt")


#Loading exonicLength info
load("DATA/exonicLength.rda")

#This is how the exonic length file was generated:
#library(GenomicFeatures)
#txdb=makeTxDbFromGFF("gencode.v19.annotation.gtf", format="gtf")
#exons.list.per.gene <- exonsBy(txdb,by="gene")
#exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
#exonicLength=as.matrix(exonic.gene.sizes )
#colnames(exonicLength)="exoniclength_bp"
#save(exonicLength, file="exonicLength.rda")


###total gene length
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#tx_by_gene <- transcriptsBy(txdb, by="gene")
#gene_lens <- max(width(tx_by_gene))



rownames(exonicLength)=transpose(strsplit(rownames(exonicLength), ".", fixed=TRUE))[[1]]

#All the genes are in the annotation file--great!
length=transpose(exonicLength)[[1]]
names(length)=rownames(exonicLength)
m=match(rownames(countData), names(length))
length=length[m]/1000
#we will use the sum of counts as library size
libsize=apply(countData,2,sum)/10^6
# TPM normalization
tpm=countData
for (j in c(1:ncol(tpm))) tpm[,j]=countData[,j]/libsize[j]
#RPKM normalization
rpkm=tpm
for (j in c(1:ncol(rpkm))) rpkm[,j]=tpm[,j]/length

write.table(rpkm, "DATA/allData_FPKM_renormalized_IV.txt", sep="\t", row.names=TRUE)
