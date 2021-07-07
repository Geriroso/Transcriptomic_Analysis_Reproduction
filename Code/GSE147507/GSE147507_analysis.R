###########################################################################
## GSE147507: R code
###########################################################################

# This is the script containing the whole code for the analysis of the GSE150316 
# dataset.

# Set up paths
wd <- "~/Desktop/TFG/Data/files/GSE147507"
GSE150wd <- "~/Desktop/TFG/Data/files/GSE150316"
setwd(wd)
dir.create(data, showWarnings = FALSE)
dir.create(results, showWarnings = FALSE)
data <- file.path(wd, "data")
results <- file.path(wd, "results")
results150 <- file.path(GSE150wd, "results")

# Load packages
library(DESeq2)
library(xlsx)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(goProfiles)
library(GenomicFeatures)

# Access the data
cts <- read.table(file = file.path(data, "GSE147507_paperdata.txt"), header = TRUE, sep = "\t", row.names = 1) 
View(cts)

###########################################################################
## NetworkAnalyst analysis
###########################################################################
# I'm going to average first the technical replicates and then sumbit them in the application
df <- data.frame(
  Case1 = cts$Series1_NHBE_SARS.CoV.2_1,
  Case2 = cts$Series1_NHBE_SARS.CoV.2_2,
  Case3 = cts$Series1_NHBE_SARS.CoV.2_3,
  NegCtrl1 = cts$Series1_NHBE_Mock_1,
  NegCtrl2 = cts$Series1_NHBE_Mock_2,
  NegCtrl3 = cts$Series1_NHBE_Mock_3)
NAME <- colnames(df)
df <- as.matrix(df)
CLASS_CvsI <- c(rep("Infected", 3), rep("Control", 3))
colnames(df) <- NULL
NetA_GSE147 <- rbind(NAME, CLASS_CvsI, df)

# As I can't see an easy way to build the needed matrix, I'm gonna add the "#" and the ":" on the .txt file
write.table(NetA_GSE147, file = file.path(results, "NetA_GSE147_averaged.txt"), sep="\t", quote=F, col.names=NA)
## NO AVERAGING by "rowMeans()" since the matrix does not work this way

###########################################################################
## DESeq2 package: NetworkAnalyst steps
###########################################################################

# Define the parameters used. Actually, the paper only specifies |FC|>2 for the
# BALF samples. I assumed the same cutoff because no other cutoff is specified.
reads_cutoff <- 4
fc_cutoff <- 2
padj_cutoff <- 0.05
variance_filter <- 15

# Access the data
cts <- read.table(file = file.path(data, "GSE147507_paperdata.txt"), header = TRUE, sep = "\t", row.names = 1)
View(cts)

df <- data.frame(
  Gid = rownames(cts),
  SARS_CoV1 = cts$Series1_NHBE_SARS.CoV.2_1,
  SARS_CoV2 = cts$Series1_NHBE_SARS.CoV.2_2,
  SARS_CoV3 = cts$Series1_NHBE_SARS.CoV.2_3,
  Mock1 = cts$Series1_NHBE_Mock_1,
  Mock2 = cts$Series1_NHBE_Mock_2,
  Mock3 = cts$Series1_NHBE_Mock_3)
rownames(df) <- rownames(cts)

# remove constant values
filter.val <- apply(df, 1, IQR, na.rm=T); #int.mat is the gene expression table
good.inx2 <- filter.val > 0;
if(sum(!good.inx2) > 0){
  df <- df[good.inx2,];
}

# Variance filter
cnt <- df[, 2:ncol(df)]
vf <- sort(apply(cnt, 1, var))
sel <- vf >= quantile(vf, p = variance_filter/100)
cnt_var <- cnt[sel,]

#Low abundance filter
af <- rowSums(cnt_var)
sel <- af > reads_cutoff
cnt_f1 <- cnt_var[sel,]
str(cnt_f1)  

#DESeq Object
all_cnt_mat <- as.matrix(cnt_f1)

condition <- factor(c(rep("SARS_CoV", 3), rep("Mock", 3)), levels = c("Mock", "SARS_CoV"))

dds <- DESeqDataSetFromMatrix(all_cnt_mat, DataFrame(condition), ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res$Gid <- rownames(res)  

res$Tag <- "NC"
res$Tag[res$log2FoldChange > log2(fc_cutoff) & res$padj < padj_cutoff] <- "Up"
res$Tag[res$log2FoldChange < -log2(fc_cutoff) & res$padj < padj_cutoff] <- "Down"
res$Tag <- factor(res$Tag, levels = c("Up", "NC", "Down"))

write.table(res, file = file.path(results, "DEA_GSE147507_NetAStep.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")


sum(res$Tag == "Up")
sum(res$Tag == "Down")  

# Metascape
## Matrix with Up&Down regulated genes in comparison with the ones obtained from the GSE150316 database
DEA_table_GSE147 <- read.table(file = file.path(results,"DEA_GSE147507_NetAStep.txt"), sep = "\t", header = TRUE, dec = ",")
DEA_table_GSE150 <- read.table(file = file.path(results150,"DEA_GSE150316_NetAStep.txt"), sep = "\t", header = TRUE, dec = ",")

DEA_table_GSE147 <- DEA_table_GSE147[, c("Gid", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Tag")]
rownames(DEA_table_GSE147) <- DEA_table_GSE147$Gid

DEA_table_GSE150 <- DEA_table_GSE150[, c("Gid", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Tag")]
rownames(DEA_table_GSE150) <- DEA_table_GSE150$Gid

up_GSE147 <- subset(DEA_table_GSE147, DEA_table_GSE147$Tag=="Up")
down_GSE147 <- subset(DEA_table_GSE147, DEA_table_GSE147$Tag=="Down")
up_GSE150 <- subset(DEA_table_GSE150, DEA_table_GSE150$Tag=="Up")
down_GSE150 <- subset(DEA_table_GSE150, DEA_table_GSE150$Tag=="Down")

# Multiple Gene List
up_genes_GSE150 <- up_GSE150$Gid
down_genes_GSE150 <- down_GSE150$Gid
up_genes_GSE147 <- up_GSE147$Gid
down_genes_GSE147 <- down_GSE147$Gid
title <- c("Up150", "Down150", "Up147", "Down147")

up_down_GSE <- matrix(0, ncol = 2, nrow = 4)
up_down_GSE <- data.frame(up_down_GSE)
up_down_GSE[,1] <- title
colnames(up_down_GSE) <- c("Names", "Genes")
up_genes_150 <- paste(up_genes_GSE150, collapse = ", ")
up_down_GSE[1,2] <- up_genes_150
down_genes_150 <- paste(down_genes_GSE150, collapse = ", ")
up_down_GSE[2,2] <- down_genes_150
up_genes_147 <- paste(up_genes_GSE147, collapse = ", ")
up_down_GSE[3,2] <- up_genes_147
down_genes_147 <- paste(down_genes_GSE147, collapse = ", ")
up_down_GSE[4,2] <- down_genes_147
View(up_down_GSE)

write.table(up_down_GSE, file = file.path(results, "Up&Down-regulated_GSE150316_and_GSE147_NetAStep.txt") ,sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



###########################################################################
## DESeq2 analysis: Default settings from the DESeq2 Package
###########################################################################

# Creation of the DESeq object
## We have already created our counts matrix in the variable "cts"
load.data <- FALSE
if(load.data){
  load(file = file.path(data, "Counts.Rda"))
  load(file = file.path(data, "Top.Rda"))
}else{
  condition <- factor(c(rep("Mock", 3), rep("SARS_CoV", 3)), levels = c("Mock", "SARS_CoV"))
  coldata <- data.frame(condition)
  rownames(coldata) <- colnames(cts)
  
  # Creation of the DESeq object
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
  
  # Pre-filtering 
  ## Here we delete the genes which do not have, at least, 10 reads in reference 
  ## to the sum of all columns
  keep <- rowSums(counts(dds)) >= 10 
  dds <- dds[keep,]
  
  # Differential Expression Analysis (DEA)
  ## The package gives indications in order to apply a Shrinkage or modifications
  ## in the analysis by modifying the p-value. I will do no changes.
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  
  # Here we have obtained the Top Table from the DESeq2
  
  # Now we will apply the cutoffs in order to classify the genes between Up&Down-regulated
  res$Tag <- "NC"
  res$Tag[res$log2FoldChange > log2(fc_cutoff) & res$padj < padj_cutoff] <- "Up"
  res$Tag[res$log2FoldChange < -log2(fc_cutoff) & res$padj < padj_cutoff] <- "Down"
  res$Tag <- factor(res$Tag, levels = c("Up", "NC", "Down"))
  
  save(dds, file = file.path(data, "Counts_DS.Rda"))
  save(res, file = file.path(data, "Top_DS.Rda"))
}

sum(res$Tag == "Up")
sum(res$Tag == "Down")

# Finally, we create our Top Table file
write.xlsx(res, file = file.path(results, "DEA_GSE147507_DS.xlsx"), sheetName = "DESeq2 data")
write.table(res, file = file.path(results, "DEA_GSE147507_DS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")

# Metascape
## Matrix with Up&Down regulated genes in comparison with the ones obtained from the GSE150316 database
DEA_table_GSE147 <- read.table(file = file.path(results,"DEA_GSE147507_copyBALF.txt"), sep = "\t", header = TRUE, dec = ",")
DEA_table_GSE150 <- read.table(file = file.path(results150,"DEA_GSE150316_copyBALF.txt"), sep = "\t", header = TRUE, dec = ",")

DEA_table_GSE147 <- DEA_table_GSE147[, c("Gid", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Tag")]
rownames(DEA_table_GSE147) <- DEA_table_GSE147$Gid

DEA_table_GSE150 <- DEA_table_GSE150[, c("Gid", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Tag")]
rownames(DEA_table_GSE150) <- DEA_table_GSE150$Gid

up_GSE147 <- subset(DEA_table_GSE147, DEA_table_GSE147$Tag=="Up")
down_GSE147 <- subset(DEA_table_GSE147, DEA_table_GSE147$Tag=="Down")
up_GSE150 <- subset(DEA_table_GSE150, DEA_table_GSE150$Tag=="Up")
down_GSE150 <- subset(DEA_table_GSE150, DEA_table_GSE150$Tag=="Down")

# Multiple Gene List
up_genes_GSE150 <- up_GSE150$Gid
down_genes_GSE150 <- down_GSE150$Gid
up_genes_GSE147 <- up_GSE147$Gid
down_genes_GSE147 <- down_GSE147$Gid
title <- c("Up150", "Down150", "Up147", "Down147")

up_down_GSE <- matrix(0, ncol = 2, nrow = 4)
up_down_GSE <- data.frame(up_down_GSE)
up_down_GSE[,1] <- title
colnames(up_down_GSE) <- c("Names", "Genes")
up_genes_150 <- paste(up_genes_GSE150, collapse = ", ")
up_down_GSE[1,2] <- up_genes_150
down_genes_150 <- paste(down_genes_GSE150, collapse = ", ")
up_down_GSE[2,2] <- down_genes_150
up_genes_147 <- paste(up_genes_GSE147, collapse = ", ")
up_down_GSE[3,2] <- up_genes_147
down_genes_147 <- paste(down_genes_GSE147, collapse = ", ")
up_down_GSE[4,2] <- down_genes_147
View(up_down_GSE)

write.table(up_down_GSE, file = file.path(results, "Up&Down-regulated_GSE150316_and_GSE147_copyBALF.txt") ,sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

###########################################################################
## Extra
###########################################################################

###########################################################################
## Excel
###########################################################################
Symbols <- rownames(cts)

#columns(org.Hs.eg.db) # returns list of available keytypes
EntrezID = mapIds(org.Hs.eg.db,
                  keys=Symbols, #Column containing Ensembl gene ids
                  column="ENTREZID",
                  keytype="SYMBOL",
                  multiVals="first")
EnsemblID = mapIds(org.Hs.eg.db,
                   keys=Symbols,
                   column="ENSEMBL",
                   keytype="SYMBOL",
                   multiVals="first")

final_mat <- cbind(Symbols, EntrezID, EnsemblID, cts)
final_mat <- merge(final_mat, res, by = 0, all = FALSE)
final_mat <- final_mat[, -1]
write.table(final_mat, file = file.path(results, "GSE147507_results_matrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")

View(final_mat)

###########################################################################
## GoProfiles
###########################################################################
up_GSE147 <- subset(res, res$Tag=="Up")
down_GSE147 <- subset(res, res$Tag=="Down")


up_GSE147$EntrezID <- mapIds(org.Hs.eg.db,
                             keys=up_GSE147$Gid, #Column containing Ensembl gene ids
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first")
up_GSE147_Entrez <- up_GSE147$EntrezID
up_GSE147_Entrez <- up_GSE147_Entrez[!is.na(up_GSE147_Entrez)]

down_GSE147$EntrezID <- mapIds(org.Hs.eg.db,
                               keys=down_GSE147$Gid, #Column containing Ensembl gene ids
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")
down_GSE147_Entrez <- down_GSE147$EntrezID
down_GSE147_Entrez <- down_GSE147_Entrez[!is.na(down_GSE147_Entrez)]

GO_Up_GSE147 <- basicProfile (up_GSE147_Entrez, onto="BP", level=2, orgPackage="org.Hs.eg.db") 
printProfiles(GO_Up_GSE147, percentage = T)
plotProfiles (GO_Up_GSE147, aTitle="Up-regulated genes GSE147")

GO_Down_GSE147 <- basicProfile (down_GSE147_Entrez, onto="BP", level=2, orgPackage="org.Hs.eg.db") 
printProfiles(GO_Down_GSE147, percentage = T)
plotProfiles (GO_Down_GSE147, percentage = T, aTitle="Down-regulated genes GSE147")

GO_Up_Down_GSE147507 <-mergeProfilesLists(GO_Up_GSE147, GO_Down_GSE147, profNames=c("Up", "Down"))
plotProfiles (GO_Up_Down_GSE147507, percentage=T,aTitle="Up vs Down", legendText = T)

fisherGOProfiles(GO_Up_GSE147$BP, GO_Down_GSE147$BP, method="BH")








