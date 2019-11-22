####Plot for Sfari&ID genes and variants enrichment per module####
library(dplyr)
library(magrittr)
library(readr)
library(gplots)
library(tidyr)
library(WGCNA)
library(vcd)

setwd ("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/artigo_expressao/ASDiPSCTranscriptome")

#pvalues and odds ratio table
modall=read.delim("DATA/Sfari_ID_Variants_Module_enrichment_matrix.txt")


#with all significant values
#getting only the numerical values from the table
modpvnum=modall[c(1:7),c(2,11,13,14,26,30)]
modornum=modall[c(1:7),c(33,42,44,45,57,61)] 
colnames(modornum)=colnames(modpvnum)


ySymbols=as.vector(modall[c(1:7),1])


#color pallete
pallete=colorRampPalette(c("white", "red"))
#vector of colors for x and y axes
colorx=c("MEturquoise", "MEblue","MEturquoise","MEtan", "MEred", "MEpurple")


####plotting table without grids (Fig.5)####
sizeGrWindow(20,12)
pdf("RESULTS/sfari_variants_enrichment_figure5.pdf")
#par(mar = c(bottom, left, up, right))
par(mar = c(10, 6, 10, 1));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = modpvnum,#pvalues
               xLabels = colorx,#colors for x-axis labeling
               xSymbols = c("MNpc1-turquoise", "MNpc10-blue", "MNeu1-turquoise", "MNeu2-tan", "MNeu14-red","MNeu18-purple"), #names for x-axis labeling
               yLabels= ySymbols, #names for y-axis labeling
               colorLabels = FALSE,
               xLabelsAdj = 0.5, #center the text label of x-axis
               xLabelsAngle = 0,
               colors = pallete(50),
               textMatrix = modornum, #paste the odds-ratio values in the table
               setStdMargins = FALSE,
               cex.text = 1, #size of pasted text in the matrix
               cex.lab.x = 0.7,
               cex.lab.y = 1,
               x.adj.lab.y = 0.5,#center the text label of y axis
               zlim = c(0,5), #set the color-coded range
               main = paste("Variants and Sfari genes enrichment per module"))
legend2 = grid_legend(0.9, 0.9,labels = "-log10(padj-value)", draw = FALSE, frame = FALSE)

grid.draw(grobTree(legend2, vp = viewport(x = 0.98, angle = 90)))
dev.off()


