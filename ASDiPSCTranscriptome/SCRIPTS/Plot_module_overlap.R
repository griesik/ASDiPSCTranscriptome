####Plotting results for Modules Overlap####
library(dplyr)
library(WGCNA)
library(magrittr)
library(readr)
library(gplots)
library(tidyr)
library(vcd)



setwd ("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/artigo_expressao/ASDiPSCTranscriptome")

#pvalues table
modpv=read.delim("DATA/module_overlap_final_matrix_pvalue.txt")
#odds-ratio table
modor=read.delim("DATA/module_overlap_final_matrix_OR.txt")
#getting only the numerical values from the tables
modpvnum=modpv[,c(3:5)]
modornum=modor[,c(3:5)]

#color pallete
pallete=colorRampPalette(c("pink", "red"))
#vector of colors for x and y axes
colorx=c("MEblue","MEturquoise","MEpurple")
colory=as.vector(modor[,7])


###ploting table without grids####
sizeGrWindow(10,6)

pdf("RESULTS/module_overlap_final_figure_nogrids2.pdf")

#par(mar = c(bottom, left, up, right))
par(mar = c(4, 8.5, 2, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = modpvnum,#pvalues
               xLabels = colorx,#colors for x-axis labeling
               xSymbols = c("MNPC10-blue", "MNeur1-turquoise", "MNeur18-purple"), #names for x-axis labeling
               yLabels = colory, ##colors for y-axis labeling
               ySymbols= modor[,2], #names for y-axis labeling
               colorLabels = TRUE,
               xLabelsAngle = 0,#set the x-axis to the horizontal position
               xLabelsAdj = 0.5, #center the text label of x-axis
               colors = pallete(150),
               naColor = "white", #NA characters should be white
               textMatrix = modornum, #paste the odds-ratio values in the table
               setStdMargins = FALSE,
               cex.text = 1, #size of pasted text in the matrix
               cex.lab.x = 0.8,
               cex.lab.y = 0.8,
               x.adj.lab.y = 0.5,#center the text label of y axis
               zlim = c(0,80), #set the color-coded range
               main = paste("Module Overlap"))

legend(x = as.numeric(0.8),y = as.numeric(1),
       bty = "n",
       legend = unique(modpv$Cell.source),
       col = c("green", "greenyellow", "yellow", "red"), 
       lty= 1,             
       lwd = 5,           
       cex=.7)


legend2 = grid_legend(0.9, 0.9,labels = "-log(padj-value)", draw = FALSE, frame = FALSE)

grid.draw(grobTree(legend2, vp = viewport(x = 0.93, angle = 90)))

dev.off()