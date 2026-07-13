#this script performs deseq2 analysis of sequence enrichment on NKX2.1 binding using seq count .txt files
#Freyda Lenihan-Geels, 2025

library(DESeq2)
library(apeglm) 
library(pheatmap)

#import.txt files of seq counts for each sample

#rename columns
colnames(C11) = c("sequence", "C11")
colnames(C12) = c("sequence", "C12")
colnames(C13) = c("sequence", "C13")
colnames(C21) = c("sequence", "C21")
colnames(C22) = c("sequence", "C22")
colnames(C23) = c("sequence", "C23")
colnames(CD1) = c("sequence", "CD1")
colnames(CD2) = c("sequence", "CD2")
colnames(CD3) = c("sequence", "CD3")

colnames(F21) = c("sequence", "F21")
colnames(F22) = c("sequence", "F22")
colnames(F23) = c("sequence", "F23")
colnames(F11) = c("sequence", "F11")
colnames(F12) = c("sequence", "F12")
colnames(F13) = c("sequence", "F13")
colnames(FD1) = c("sequence", "FD1")
colnames(FD2) = c("sequence", "FD2")
colnames(FD3) = c("sequence", "FD3")

colnames(A51) = c("sequence", "A51")
colnames(A52) = c("sequence", "A52")
colnames(A53) = c("sequence", "A53")
colnames(A21) = c("sequence", "A21")
colnames(A22) = c("sequence", "A22")
colnames(A23) = c("sequence", "A23")
colnames(A11) = c("sequence", "A11")
colnames(A12) = c("sequence", "A12")
colnames(A13) = c("sequence", "A13")
colnames(AD1) = c("sequence", "AD1")
colnames(AD2) = c("sequence", "AD2")
colnames(AD3) = c("sequence", "AD3")

#merge dataframes
dat_flank = merge(FD1, FD2, by = "sequence", all = TRUE)
dat_flank = merge(dat, FD3, by = "sequence", all = TRUE)
dat_flank = merge(dat, F21, by = "sequence", all = TRUE)
dat_flank = merge(dat, F22, by = "sequence", all = TRUE)
dat_flank = merge(dat, F23, by = "sequence", all = TRUE)
dat_flank = merge(dat, F11, by = "sequence", all = TRUE)
dat_flank = merge(dat, F12, by = "sequence", all = TRUE)
dat_flank = merge(dat, F13, by = "sequence", all = TRUE)

dat_core = merge(CD1, CD2, by = "sequence", all = TRUE)
dat_core = merge(dat_core, CD3, by = "sequence", all = TRUE)
dat_core = merge(dat_core, C21, by = "sequence", all = TRUE)
dat_core = merge(dat_core, C22, by = "sequence", all = TRUE)
dat_core = merge(dat_core, C23, by = "sequence", all = TRUE)
dat_core = merge(dat_core, C11, by = "sequence", all = TRUE)
dat_core = merge(dat_core, C12, by = "sequence", all = TRUE)
dat_core = merge(dat_core, C13, by = "sequence", all = TRUE)

dat_all = merge(AD1, AD2, by = "sequence", all = TRUE)
dat_all = merge(dat_all, AD3, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A11, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A12, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A13, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A21, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A22, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A23, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A51, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A52, by = "sequence", all = TRUE)
dat_all = merge(dat_all, A53, by = "sequence", all = TRUE)

#FLANK LIBRARY ONLY - remove wildtype sequence due to contamination (of CORE library seqs)
dat_flank <- subset(dat_flank, sequence != "CCAGTCAAGTGTTC")

rownames(dat) <- make.unique(as.character(dat$sequence))
dat <- dat[,-1]
colnames(dat) = c("Unbound_rep1", "Unbound_rep2", "Unbound_rep3", "0.2_rep1", "0.2_rep2", "0.2_rep3", "0.1_rep1", "0.1_rep2", "0.1_rep3")
dat[is.na(dat)] = 0
dat <- dat[dat$Unbound_rep1 >= 5 & dat$Unbound_rep2 >= 5 & dat$Unbound_rep3 >= 5, ] #removes any sequences with under 5 counts in unbound libraries

# Check if the row names in coldata match the column names in dat
all(rownames(coldata) == colnames(dat))

#Make sure conditions are set for unbound vs the bound samples
coldata$Condition <- factor(coldata$Condition, levels = c("Unbound", "0.2", "0.1")) # add 0.5 for all library

# Create a DESeqDataSet object for differential expression analysis and perform analysis
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design = ~ Condition)
dds<- DESeq(dds)

# Get the differential expression results
res <- results(dds)

# Get a specific result
res_0.1 = results(dds, name="Condition_0.1_vs_Unbound")
res_0.2 = results(dds, name="Condition_0.2_vs_Unbound")
res_0.5 = results(dds, name0"Condition_0.5_vs_Unbound") #all library only

# Shrink the log-fold changes using apeglm method
res_0.1 <- lfcShrink(dds, coef = "Condition_0.1_vs_Unbound", type = "apeglm")
res_0.2 <- lfcShrink(dds, coef = "Condition_0.2_vs_Unbound", type = "apeglm")
res_0.5 <- lfcShrink(dds, coef = "Condition_0.5_vs_Unbound", type = "apeglm") #all library only

# Plot the MA plot of shrunken results
#png("MAplot1.png", res = 600, width = 4500, height = 2000)
plotMA(
  res_0.1,
  colNonSig = "gray",
  colSig = "darkorchid1",
  ylim = c(-6, 2),
  cex = 0.8,
  cex.axis = 1.2,
  cex.lab = 1.3,
  cex.main = 1.4
)
dev.off()

png("MAplot2.png", res = 600, width = 4500, height = 2000)
plotMA(
  res_0.2,
  colNonSig = "gray",
  colSig = "darkorchid1",
  ylim = c(-1, 5),
  cex = 0.8,
  cex.axis = 1.2,
  cex.lab = 1.3,
  cex.main = 1.4
)
dev.off()

vsd <- varianceStabilizingTransformation(dds, blind=FALSE) #for core library
rld <- rlog(dds, blind = FALSE) #for all and flank libraries

plotPCA(vsd, intgroup=c("Condition","Replicate"))
plotPCA(rld, intgroup=c("Condition","Replicate"))

# Order the results by log-fold change in decreasing order
resOrdered_0.1 <- res_0.1[order(res_0.1$log2FoldChange, decreasing = TRUE),]

# Write the ordered results to a CSV file
write.csv(as.data.frame(resOrdered_0.1), file = "/core_0.1_results.csv")

#Create Euclidean sample-sample distance plots. Use vsd or rlg according to library
# Calculate sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

# Convert to matrix for heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

# Choose a color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Plot the heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Distance Heatmap")

