#####ID and Sfari genes enrichment ####
#####NPC enrichment per module####
#select only protein coding genes for this analysis
codinggene=subset(NPCresult, NPCresult$Gene_type=="protein_coding")


#calculate the total number of genes
tg=as.numeric(nrow(codinggene)) #total number of genes
idg=codinggene[!is.na(codinggene$ID_score),]
idg=as.numeric(nrow(idg)) #total number of ID genes expressed in the data set
idog=codinggene[!is.na(codinggene$ID_score)&is.na(codinggene$Sfari_score),]
idog=as.numeric(nrow(idog)) #total number of ID exclusive (only) genes expressed in the data set
sfg=codinggene[!is.na(codinggene$Sfari_score),]
sfg=as.numeric(nrow(sfg)) #total number of sfari genes expressed in the dataset
sfog=codinggene[!is.na(codinggene$Sfari_score)&is.na(codinggene$ID_score),]
sfog=as.numeric(nrow(sfog)) #total number of sfari exclusive (only) genes expressed in the dataset
idsfg=codinggene[!is.na(codinggene$Sfari_score)&!is.na(codinggene$ID_score),]
idsfg=as.numeric(nrow(idsfg)) #total number of sfari&ID (overlap) genes expressed in the dataset


#create the table for fisher test for ID all
enrichtabid<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                        nmnsg=numeric(0), module_genes=numeric(0),
                        total_genes=numeric(0),pvalue=numeric(0), 
                        oddsratio=numeric(0),
                        conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$ID_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=idg-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabid[j,c(1:4)]=l
enrichtabid[j,5]=as.numeric(mg)
enrichtabid[j,6]=as.numeric(tg)}

#calculate p-values
for(j in c(1:nrow(enrichtabid)))
{counts = (matrix(data = as.numeric(as.character(enrichtabid[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabid$pvalue[j]=res$p.value
enrichtabid$oddsratio[j]=res$estimate
enrichtabid$conf_intervalmin[j]=res$conf.int[1]
enrichtabid$conf_intervalmax[j]=res$conf.int[2]
}


#create the table for fisher test for ID only
enrichtabido<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                         nmnsg=numeric(0), module_genes=numeric(0),
                         total_genes=numeric(0),pvalue=numeric(0), 
                         oddsratio=numeric(0),
                         conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$ID_score)&is.na(modj$Sfari_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=idog-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabido[j,c(1:4)]=l
enrichtabido[j,5]=as.numeric(mg)
enrichtabido[j,6]=as.numeric(tg)}

#calculate p-values
for(j in c(1:nrow(enrichtabido)))
{counts = (matrix(data = as.numeric(as.character(enrichtabido[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabido$pvalue[j]=res$p.value
enrichtabido$oddsratio[j]=res$estimate
enrichtabido$conf_intervalmin[j]=res$conf.int[1]
enrichtabido$conf_intervalmax[j]=res$conf.int[2]
}



#create the table for fisher test for Sfari genes all
enrichtabsf<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                        nmnsg=numeric(0), module_genes=numeric(0),
                        total_genes=numeric(0),pvalue=numeric(0), 
                        oddsratio=numeric(0),
                        conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$Sfari_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=sfg-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabsf[j,c(1:4)]=l
enrichtabsf[j,5]=as.numeric(mg)
enrichtabsf[j,6]=as.numeric(tg)}


for(j in c(1:nrow(enrichtabsf)))
{counts = (matrix(data = as.numeric(as.character(enrichtabsf[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabsf$pvalue[j]=res$p.value
enrichtabsf$oddsratio[j]=res$estimate
enrichtabsf$conf_intervalmin[j]=res$conf.int[1]
enrichtabsf$conf_intervalmax[j]=res$conf.int[2]
}


#create the table for fisher test for Sfari genes only
enrichtabsfo<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                         nmnsg=numeric(0), module_genes=numeric(0),
                         total_genes=numeric(0),pvalue=numeric(0), 
                         oddsratio=numeric(0),
                         conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$Sfari_score)&is.na(modj$ID_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=sfog-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabsfo[j,c(1:4)]=l
enrichtabsfo[j,5]=as.numeric(mg)
enrichtabsfo[j,6]=as.numeric(tg)}


for(j in c(1:nrow(enrichtabsfo)))
{counts = (matrix(data = as.numeric(as.character(enrichtabsfo[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabsfo$pvalue[j]=res$p.value
enrichtabsfo$oddsratio[j]=res$estimate
enrichtabsfo$conf_intervalmin[j]=res$conf.int[1]
enrichtabsfo$conf_intervalmax[j]=res$conf.int[2]
}


#create the table for fisher test for Sfari & ID overlap genes
enrichtabidsf<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                          nmnsg=numeric(0), module_genes=numeric(0),
                          total_genes=numeric(0),pvalue=numeric(0), 
                          oddsratio=numeric(0),
                          conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$Sfari_score)&!is.na(modj$ID_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=idsfg-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabidsf[j,c(1:4)]=l
enrichtabidsf[j,5]=as.numeric(mg)
enrichtabidsf[j,6]=as.numeric(tg)}


for(j in c(1:nrow(enrichtabidsf)))
{counts = (matrix(data = as.numeric(as.character(enrichtabidsf[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabidsf$pvalue[j]=res$p.value
enrichtabidsf$oddsratio[j]=res$estimate
enrichtabidsf$conf_intervalmin[j]=res$conf.int[1]
enrichtabidsf$conf_intervalmax[j]=res$conf.int[2]
}

#combining both tables and calculate padj value
#ID all
enrichtabid[,9]="ID_all"
enrichtabid[,10]=paste(enrichtabid[,9],rownames(enrichtabid))
rownames(enrichtabid)=enrichtabid[,10]
enrichtabid=enrichtabid[,c(1:8)]
enrichtabid=enrichtabid[c(1,2,5:12,3,4),]

#ID only
enrichtabido[,9]="ID_only"
enrichtabido[,10]=paste(enrichtabido[,9],rownames(enrichtabido))
rownames(enrichtabido)=enrichtabido[,10]
enrichtabido=enrichtabido[,c(1:8)]
enrichtabido=enrichtabido[c(1,2,5:12,3,4),]

#Sfari all
enrichtabsf[,9]="Sfari_all"
enrichtabsf[,10]=paste(enrichtabsf[,9],rownames(enrichtabsf))
rownames(enrichtabsf)=enrichtabsf[,10]
enrichtabsf=enrichtabsf[,c(1:8)]
enrichtabsf=enrichtabsf[c(1,2,5:12,3,4),]

#Sfari_only
enrichtabsfo[,9]="Sfari_only"
enrichtabsfo[,10]=paste(enrichtabsfo[,9],rownames(enrichtabsfo))
rownames(enrichtabsfo)=enrichtabsfo[,10]
enrichtabsfo=enrichtabsfo[,c(1:8)]
enrichtabsfo=enrichtabsfo[c(1,2,5:12,3,4),]

#Sfari&ID overlap
enrichtabidsf[,9]="Sfari&ID"
enrichtabidsf[,10]=paste(enrichtabidsf[,9],rownames(enrichtabidsf))
rownames(enrichtabidsf)=enrichtabidsf[,10]
enrichtabidsf=enrichtabidsf[,c(1:8)]
enrichtabidsf=enrichtabidsf[c(1,2,5:12,3,4),]

#combine all
enrichtabNPC=rbind(enrichtabsf,enrichtabsfo,enrichtabid,enrichtabido,enrichtabidsf)
enrichtabNPC[,9]="NPC"
enrichtabNPC[,10]=paste(enrichtabNPC[,9],rownames(enrichtabNPC))
rownames(enrichtabNPC)=enrichtabNPC[,10]
enrichtabNPC=enrichtabNPC[,c(1:8)]




#####Neuron enrichment per module####
#select only protein coding genes for this analysis
Neuronresult=read.csv("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/Neurons/dreamResults_Neuron_nDEGs_HKgenes_p0.4_k4.csv")
codinggene=subset(Neuronresult, Neuronresult$Gene_type=="protein_coding")


#calculate the total number of genes
tg=as.numeric(nrow(codinggene)) #total number of genes
idg=codinggene[!is.na(codinggene$ID_score),]
idg=as.numeric(nrow(idg)) #total number of ID genes expressed in the data set
idog=codinggene[!is.na(codinggene$ID_score)&is.na(codinggene$Sfari_score),]
idog=as.numeric(nrow(idog)) #total number of ID exclusive (only) genes expressed in the data set
sfg=codinggene[!is.na(codinggene$Sfari_score),]
sfg=as.numeric(nrow(sfg)) #total number of sfari genes expressed in the dataset
sfog=codinggene[!is.na(codinggene$Sfari_score)&is.na(codinggene$ID_score),]
sfog=as.numeric(nrow(sfog)) #total number of sfari exclusive (only) genes expressed in the dataset
idsfg=codinggene[!is.na(codinggene$Sfari_score)&!is.na(codinggene$ID_score),]
idsfg=as.numeric(nrow(idsfg)) #total number of sfari&ID (overlap) genes expressed in the dataset


#create the table for fisher test for ID all
enrichtabid<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                        nmnsg=numeric(0), module_genes=numeric(0),
                        total_genes=numeric(0),pvalue=numeric(0), 
                        oddsratio=numeric(0),
                        conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$ID_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=idg-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabid[j,c(1:4)]=l
enrichtabid[j,5]=as.numeric(mg)
enrichtabid[j,6]=as.numeric(tg)}

#calculate p-values
for(j in c(1:nrow(enrichtabid)))
{counts = (matrix(data = as.numeric(as.character(enrichtabid[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabid$pvalue[j]=res$p.value
enrichtabid$oddsratio[j]=res$estimate
enrichtabid$conf_intervalmin[j]=res$conf.int[1]
enrichtabid$conf_intervalmax[j]=res$conf.int[2]
}


#create the table for fisher test for ID only
enrichtabido<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                         nmnsg=numeric(0), module_genes=numeric(0),
                         total_genes=numeric(0),pvalue=numeric(0), 
                         oddsratio=numeric(0),
                         conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$ID_score)&is.na(modj$Sfari_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=idog-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabido[j,c(1:4)]=l
enrichtabido[j,5]=as.numeric(mg)
enrichtabido[j,6]=as.numeric(tg)}

#calculate p-values
for(j in c(1:nrow(enrichtabido)))
{counts = (matrix(data = as.numeric(as.character(enrichtabido[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabido$pvalue[j]=res$p.value
enrichtabido$oddsratio[j]=res$estimate
enrichtabido$conf_intervalmin[j]=res$conf.int[1]
enrichtabido$conf_intervalmax[j]=res$conf.int[2]
}



#create the table for fisher test for Sfari genes all
enrichtabsf<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                        nmnsg=numeric(0), module_genes=numeric(0),
                        total_genes=numeric(0),pvalue=numeric(0), 
                        oddsratio=numeric(0),
                        conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$Sfari_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=sfg-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabsf[j,c(1:4)]=l
enrichtabsf[j,5]=as.numeric(mg)
enrichtabsf[j,6]=as.numeric(tg)}


for(j in c(1:nrow(enrichtabsf)))
{counts = (matrix(data = as.numeric(as.character(enrichtabsf[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabsf$pvalue[j]=res$p.value
enrichtabsf$oddsratio[j]=res$estimate
enrichtabsf$conf_intervalmin[j]=res$conf.int[1]
enrichtabsf$conf_intervalmax[j]=res$conf.int[2]
}


#create the table for fisher test for Sfari genes only
enrichtabsfo<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                         nmnsg=numeric(0), module_genes=numeric(0),
                         total_genes=numeric(0),pvalue=numeric(0), 
                         oddsratio=numeric(0),
                         conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$Sfari_score)&is.na(modj$ID_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=sfog-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabsfo[j,c(1:4)]=l
enrichtabsfo[j,5]=as.numeric(mg)
enrichtabsfo[j,6]=as.numeric(tg)}


for(j in c(1:nrow(enrichtabsfo)))
{counts = (matrix(data = as.numeric(as.character(enrichtabsfo[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabsfo$pvalue[j]=res$p.value
enrichtabsfo$oddsratio[j]=res$estimate
enrichtabsfo$conf_intervalmin[j]=res$conf.int[1]
enrichtabsfo$conf_intervalmax[j]=res$conf.int[2]
}


#create the table for fisher test for Sfari & ID overlap genes
enrichtabidsf<-data.frame(msg=numeric(0),mnsg=numeric(0),nmsg=numeric(0),
                          nmnsg=numeric(0), module_genes=numeric(0),
                          total_genes=numeric(0),pvalue=numeric(0), 
                          oddsratio=numeric(0),
                          conf_intervalmin=numeric(0), conf_intervalmax=numeric(0))

#calculate and place the numbers in the table
levelsj=codinggene$moduleLabel
for (j in levels(levelsj))
{modj=subset(codinggene, codinggene$moduleLabel==j)
mg=as.numeric(nrow(modj)) #number of genes in the module
msg=modj[!is.na(modj$Sfari_score)&!is.na(modj$ID_score),]
msg=as.numeric(nrow(msg))#number of genes that are in the module and in the subset of interest (ID or Sfari)
mnsg=mg-msg#number of genes that are in the module but not in the subset of interest (ID or Sfari)
nmsg=idsfg-msg #number of genes that are NOT the module but that are in the subset of interest (ID or Sfari)
nmnsg=tg-(msg+mnsg+nmsg) #number of genes that are not in the module and not in the subset of interest (ID or Sfari)
l=t(as.matrix(c(msg,mnsg,nmsg,nmnsg)))
enrichtabidsf[j,c(1:4)]=l
enrichtabidsf[j,5]=as.numeric(mg)
enrichtabidsf[j,6]=as.numeric(tg)}


for(j in c(1:nrow(enrichtabidsf)))
{counts = (matrix(data = as.numeric(as.character(enrichtabidsf[j,c(1:4)])), nrow=2))
res = fisher.test(counts)
enrichtabidsf$pvalue[j]=res$p.value
enrichtabidsf$oddsratio[j]=res$estimate
enrichtabidsf$conf_intervalmin[j]=res$conf.int[1]
enrichtabidsf$conf_intervalmax[j]=res$conf.int[2]
}

#combining both tables and calculate padj value
#ID all
enrichtabid[,9]="ID_all"
enrichtabid[,10]=paste(enrichtabid[,9],rownames(enrichtabid))
rownames(enrichtabid)=enrichtabid[,10]
enrichtabid=enrichtabid[,c(1:8)]
enrichtabid=enrichtabid[c(1,2,13,15:21,3:12,14),]

#ID only
enrichtabido[,9]="ID_only"
enrichtabido[,10]=paste(enrichtabido[,9],rownames(enrichtabido))
rownames(enrichtabido)=enrichtabido[,10]
enrichtabido=enrichtabido[,c(1:8)]
enrichtabido=enrichtabido[c(1,2,13,15:21,3:12,14),]

#Sfari all
enrichtabsf[,9]="Sfari_all"
enrichtabsf[,10]=paste(enrichtabsf[,9],rownames(enrichtabsf))
rownames(enrichtabsf)=enrichtabsf[,10]
enrichtabsf=enrichtabsf[,c(1:8)]
enrichtabsf=enrichtabsf[c(1,2,13,15:21,3:12,14),]

#Sfari_only
enrichtabsfo[,9]="Sfari_only"
enrichtabsfo[,10]=paste(enrichtabsfo[,9],rownames(enrichtabsfo))
rownames(enrichtabsfo)=enrichtabsfo[,10]
enrichtabsfo=enrichtabsfo[,c(1:8)]
enrichtabsfo=enrichtabsfo[c(1,2,13,15:21,3:12,14),]

#Sfari&ID overlap
enrichtabidsf[,9]="Sfari&ID"
enrichtabidsf[,10]=paste(enrichtabidsf[,9],rownames(enrichtabidsf))
rownames(enrichtabidsf)=enrichtabidsf[,10]
enrichtabidsf=enrichtabidsf[,c(1:8)]
enrichtabidsf=enrichtabidsf[c(1,2,13,15:21,3:12,14),]

#combine all
enrichtabNeuron=rbind(enrichtabsf,enrichtabsfo,enrichtabid,enrichtabido,enrichtabidsf)
enrichtabNeuron[,9]="Neuron"
enrichtabNeuron[,10]=paste(enrichtabNeuron[,9],rownames(enrichtabNeuron))
rownames(enrichtabNeuron)=enrichtabNeuron[,10]
enrichtabNeuron=enrichtabNeuron[,c(1:8)]


#####combine NPC and Neuron data and calculate padj value####
enrichtab=rbind(enrichtabNPC,enrichtabNeuron)
#calculate padj value
p= enrichtab[,7]
FDR = as.data.frame(p.adjust(p, method = "BH", n = length(p)))
enrichtab = cbind (enrichtab,FDR)
colnames(enrichtab)[9]="FDR"

write.csv(enrichtab,file="C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/RESULTS/Sfari_ID_genes_enrichment.csv")

