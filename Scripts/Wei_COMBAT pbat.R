#Normalization, PCA and Screeplot
#Wei Li (mouse model) - PBAT

library("preprocessCore")
library("readxl")
library("sva")

#Parameters
cts <- read.table(file = "../../../Wei_counts.tsv", sep = '\t', header = TRUE)

#Template
codes = read_excel("../../../RNAseq Ben Constant .xlsx", sheet = 1)
codes = codes[ which(codes$Assay == "PBAT"), ]
codes <- codes[codes$Assay == "PBAT" & !codes$ID %in% c("A3_2", "A10_2", "A13_2", "B18_2") & !codes$Treatment_ %in% c("BMP10", "BMP8B"), ]
cts = cts[ , codes$ID]

table(is.element(colnames(cts), codes$ID))

# 
# unique(codes$`Tag of the Experimental Group`)
CTRL_pbat = which(codes$Treatment_ == "Control" & codes$Assay == "PBAT")
BMP4_pbat = which(codes$Treatment_ == "BMP4" & codes$Assay == "PBAT")
BMP9_pbat = which(codes$Treatment_ == "BMP9" & codes$Assay == "PBAT")

#Do the process for the pbat only
cts_pbat = cts[ , which(codes$Assay == "PBAT")]
rm(cts)
cts_pbat = data.matrix(cts_pbat)
cts_pbat <- cts_pbat[(rowSums(cts_pbat)>dim(cts_pbat)[2]),] #Exlude low expressed counts


# #CPM (Counts Per Million)
# #BiocManager::install("edgeR")
# library("edgeR")
# CPM = cpm(cts_pbat)
# write.table(CPM, sep = "\t", file ="Wei_CPMs.txt")


#Normalization
nrm=normalize.quantiles(log2(1+cts_pbat))
dimnames(nrm) = dimnames(cts_pbat)




#COMBAT batch correction
mod_ = model.matrix(~as.factor(codes$Treatment_))
corrected_counts <- ComBat(dat=as.matrix(nrm), batch=codes$BatchID, mod=mod_, par.prior=TRUE, prior.plots=FALSE)
write.csv(corrected_counts, file = "Batch Corrected Counts (PBAT).csv")
nrm = corrected_counts




#PCA Plots
pca = prcomp(t(nrm))
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #How much variation in the original data each PC accounts for

#Control and Diet points for the PCA plot
condition_ = c(1:dim(cts_pbat)[2])
condition_[CTRL_pbat]="black"
condition_[BMP4_pbat]="#e34a33"
condition_[BMP9_pbat]="#fdbb84"


pdf(file="Wei_PCA_pbat_after_batch_correction.pdf", width = 8, height = 8)

# Increase the left margin for the plot
par(mar = c(5, 6, 4, 2) + 0.1, xpd = TRUE)

plot(pca$x[,1], pca$x[,2], col=condition_, 
     xlab=paste("PC1", " - ", pca.var.per[1], "%"), 
     ylab=paste("PC2", " - ", pca.var.per[2], "%"), 
     main = "pbat", 
     pch=20, 
     cex=3.6, 
     cex.lab=2,  # Increased font size for axis labels
     cex.axis=1.8, # Increased font size for axis ticks
     cex.main=3,   # Increased font size for main title
     cex.sub=1.8,  # Increased font size for subtitle
     las=1)

legend(x="topleft", 
       cex=1.8,      # Increased font size for legend text
       pt.cex=2,     # Increased size for legend points
       legend=unique(c("CTRL", "BMP4", "BMP9")), 
       col=c("black", "#e34a33", "#fdbb84"), 
       pch=20)

#text(pca$x[,1], pca$x[,2], colnames(nrm), pos = 1, offset=.3, cex=.4)
dev.off()

pdf(file="Wei_ScreePlot pbat.pdf")
barplot(pca.var.per, main="Scree Plot", xlab="Principal Components", ylab="Percent Variation", names.arg = pca.var.per)
dev.off()


