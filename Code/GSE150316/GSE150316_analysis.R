###########################################################################
## GSE150316: R code
###########################################################################

# This is the script containing the whole code for the analysis of the GSE150316 
# dataset.

# Set up paths
wd <- "~/Desktop/TFG/Data/files/GSE150316"
setwd(wd)
dir.create(dat, showWarnings = FALSE)
dir.create(res, showWarnings = FALSE)
data <- file.path(wd, "data")
results <- file.path(wd, "results")

# Load packages
library(DESeq2)
library(xlsx)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(goProfiles)
library(GenomicFeatures)

# Access the data
cts <- read.table(file = file.path(data, "GSE150316_paperdata.txt"), header = TRUE, sep = "\t", row.names = 1) 

###########################################################################
## NetworkAnalyst analysis
###########################################################################
# I'm going to average first the technical replicates and then sumbit them in the application
df <- data.frame(
  Case1 = rowSums(cts[, c("case1.lung1", "case1.lung2", "case1.lung3", "case1.lung4")]),
  Case2 = rowSums(cts[, c("case2.lung1", "case2.lung2", "case2.lung3")]),
  Case3 = rowSums(cts[, c("case3.lung1", "case3.lung2")]),
  Case4 = rowSums(cts[, c("case4.lung1", "case4.lung2")]),
  Case5 = rowSums(cts[, c("case5.lung1", "case5.lung2", "case5.lung3", "case5.lung4", "case5.lung5")]),
  NegCtrl1 = cts$NegControl1,
  NegCtrl2 = cts$NegControl2,
  NegCtrl3 = cts$NegControl3,
  NegCtrl4 = cts$NegControl4,
  NegCtrl5 = cts$NegControl5)
NAME <- colnames(df)
df <- as.matrix(df)
CLASS_CvsI <- c(rep("Infected", 5), rep("Control", 5))
colnames(df) <- NULL
NetA_GSE150 <- rbind(NAME, CLASS_CvsI, df)

# As I can't see an easy way to build the needed matrix, I'm gonna add the "#" and the ":" on the .txt file
write.table(NetA_GSE150, file = file.path(results, "NetA_GSE150_averaged.txt"), sep="\t", quote=F, col.names=NA)
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

df <- data.frame(
  Symbol = rownames(cts),
  Case1 = rowSums(cts[, c("case1.lung1", "case1.lung2", "case1.lung3", "case1.lung4")]),
  Case2 = rowSums(cts[, c("case2.lung1", "case2.lung2", "case2.lung3")]),
  Case3 = rowSums(cts[, c("case3.lung1", "case3.lung2")]),
  Case4 = rowSums(cts[, c("case4.lung1", "case4.lung2")]),
  Case5 = rowSums(cts[, c("case5.lung1", "case5.lung2", "case5.lung3", "case5.lung4", "case5.lung5")]),
  NegCtrl1 = cts$NegControl1,
  NegCtrl2 = cts$NegControl2,
  NegCtrl3 = cts$NegControl3,
  NegCtrl4 = cts$NegControl4,
  NegCtrl5 = cts$NegControl5)

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

condition <- factor(c(rep("Case", 5), rep("Ctrl", 5)), levels = c("Ctrl", "Case"))

dds <- DESeqDataSetFromMatrix(all_cnt_mat, DataFrame(condition), ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res$Gid <- rownames(res)  

res$Tag <- "NC"
res$Tag[res$log2FoldChange > log2(fc_cutoff) & res$padj < padj_cutoff] <- "Up"
res$Tag[res$log2FoldChange < -log2(fc_cutoff) & res$padj < padj_cutoff] <- "Down"
res$Tag <- factor(res$Tag, levels = c("Up", "NC", "Down"))

write.table(res, file = file.path(results, "DEA_GSE150316_NetAStep.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")

sum(res$Tag == "Up")
sum(res$Tag == "Down")  

#### Metascape Analysis
## We need to apply the requiered format from Metascape
up_down <- matrix(0, ncol = 2, nrow = 2)
up_down <- data.frame(up_down)
up_down[,1] <- title
colnames(up_down) <- c("Names", "Genes")
up_genes <- paste(up_GSE150$Gid, collapse = ", ")
up_down[1,2] <- up_genes
down_genes <- paste(down_GSE150$Gid, collapse = ", ")
up_down[2,2] <- down_genes
write.table(up_down, file = file.path(results,"Up&Down-regulated_NetAStep_GSE150316.txt"),sep = "\t", quote = F, row.names = F, col.names = T)

###########################################################################
## GSE150316 DESeq2 analysis: Default settings from the DESeq2 Package
###########################################################################

# Creation of the DESeq object
## The first time I used this package I created an Excel containing the data
## from the columns as I needed to have an idea of what I was comparing and
## what I was collapsing
## Here I open my file and then I arrange it so that it has the appropiate 
## format for the DESeq2 analysis
load.data <- FALSE
if(load.data){
  load(file = file.path(data, "Counts.Rda"))
  load(file = file.path(data, "Top.Rda"))
}else{
  coldata <- read.table(file = file.path(data, "Coldata_GSE150316_paperdata.txt"), header = T, sep = "\t", row.names = 1)
  rownames(coldata) <- gsub("-", ".", rownames(coldata))
  coldata <- coldata[,c("condition","type", "replicate")]
  coldata$condition <- factor(coldata$condition, levels = c("control", "case"))
  coldata$type <- factor(coldata$type)
  coldata$replicate <- factor(coldata$replicate, levels = c(paste0("sample",1:10)))
  
  all(rownames(coldata) %in% colnames(cts))
  all(rownames(coldata) == colnames(cts))
  
  # Creation of the DESeq object
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
  
  # Pre-filtering 
  ## Here we delete the genes which do not have, at least, 10 reads in reference 
  ## to the sum of all columns
  keep <- rowSums(counts(dds)) >= 10 
  dds <- dds[keep,]
  
  # Collapse the technical replicates
  ## This is the tool provided by the package which is specially indicated to
  ## collapse the technical replicates. In fact, it just sums up read counts
  ## based on the indications from the data of the columns, so it does the 
  ## same job as the "rowSums()" function.
  dds$run <- paste0("run", 1:21)
  ddsColl <- collapseReplicates(dds, dds$replicate, dds$run, renameCols = F)
  colnames(ddsColl[,6:10]) <- paste0("Case",1:5) 
  colData(ddsColl)
  colnames(ddsColl)
  
  # Differential Expression Analysis (DEA)
  ## Here, I conduct the DEA with the collapsed DESeq object
  ## The package gives indications in order to apply a Shrinkage or modifications
  ## in the analysis by modifying the p-value. I will do no changes.
  dds <- DESeq(ddsColl)
  res <- results(dds)
  
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
write.xlsx(res, file = file.path(results, "DEA_GSE150316_DS.xlsx"), sheetName = "DESeq2 data")
write.table(res, file = file.path(results, "DEA_GSE150316_DS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")

# Metascape Analysis
## We need to prepare the DEGs in order to introduce them on the Metascape software
## Here we keep the genes that have been tagged depending on their FC and p-adj
up_GSE150 <- subset(res, res$Tag=="Up")
down_GSE150 <- subset(res, res$Tag=="Down")
title <- c("Up", "Down")

up_down <- matrix(0, ncol = 2, nrow = 2)
up_down <- data.frame(up_down)
up_down[,1] <- title
colnames(up_down) <- c("Names", "Genes")
up_genes <- paste(rownames(up_GSE150), collapse = ", ")
up_down[1,2] <- up_genes
down_genes <- paste(rownames(down_GSE150), collapse = ", ")
up_down[2,2] <- down_genes
View(up_down)
write.table(up_down, file = file.path(results,"Up&Down-regulated_DS_GSE150316.txt"),sep = "\t", quote = F, row.names = F, col.names = T)

###########################################################################
## Extra
###########################################################################

# TPM
# Extract gene length


txdb <- makeTxDbFromGFF(file = file.path(data, "gencode.v38.annotation.gtf"),format="gtf")
txdb <- makeTxDbFromGFF(file = file.path(data, "GCF_009858895.2_ASM985889v3_genomic.gff"),format="gff")
exons.list.per.gene <- exonsBy(txdb,by="gene")
length(exons.list.per.gene)
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
exonic.gene.sizes <- as.data.frame(exonic.gene.sizes)
colnames(exonic.gene.sizes) <- "Gene_length"
exonic.gene.sizes$Gid <- rownames(exonic.gene.sizes)
View(exonic.gene.sizes)
exonic.gene.sizes <- exonic.gene.sizes[, -2]

genes <- merge(gene_info, exonic.gene.sizes, all=FALSE)
View(genes)
genes <- genes[,-1]
colnames(genes) <- c("Gid","Gene_length")

# Aquest TPM està molt mal fet, ja que no s'ha utilitzat el mateix GTF que han utilitzat al paper, s'ha utilitzat un GTF
# del SARS-CoV-2 del NCBI (a la finestra d'inici d'info pel COVID19 està allà, tota la info del SARS-CoV-2)

cts <- read.table(file = file.path(dat, "GSE150316_paperdata.txt"), header = TRUE, sep = "\t", row.names = 1) 
df <- data.frame(
  Gid = rownames(cts),
  Case1 = rowSums(cts[, c("case1.lung1", "case1.lung2", "case1.lung3", "case1.lung4")]),
  Case2 = rowSums(cts[, c("case2.lung1", "case2.lung2", "case2.lung3")]),
  Case3 = rowSums(cts[, c("case3.lung1", "case3.lung2")]),
  Case4 = rowSums(cts[, c("case4.lung1", "case4.lung2")]),
  Case5 = rowSums(cts[, c("case5.lung1", "case5.lung2", "case5.lung3", "case5.lung4", "case5.lung5")]),
  NegCtrl1 = cts$NegControl1,
  NegCtrl2 = cts$NegControl2,
  NegCtrl3 = cts$NegControl3,
  NegCtrl4 = cts$NegControl4,
  NegCtrl5 = cts$NegControl5)
View(df_cnt)
df <- merge(df, genes, all = FALSE)
df<- df[!duplicated(df$Gid),]

df_cnt <- df[, c("Case1", "Case2", "Case3", "Case4", "Case5", "NegCtrl1", "NegCtrl2", "NegCtrl3", "NegCtrl4", "NegCtrl5")]
rownames(df_cnt) <- df$Gid

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## En principi ho és, però seria necessari assegurar-se que l'ordre dels gens és el mateix per no barrejar gene lengths

tpms <- apply(df_cnt, 2, function(x) tpm(x, df$Gene_length))
View(tpms)
write.table(tpms, file = file.path(results, "TPM_counts_GSE150316.txt"), sep = "\t", quote = F, col.names = NA)

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
final_mat <-merge(final_mat, expr_res, by = 0, all = FALSE)
final_mat <- final_mat[, -1]
write.table(final_mat, file = file.path(results, "GSE150316_results_matrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")

###########################################################################
## GoProfiles
###########################################################################
# This was done with the DESeq2 package: NetworkAnalyst steps
up_GSE150 <- subset(res, res$Tag=="Up")
down_GSE150 <- subset(res, res$Tag=="Down")


up_GSE150$EntrezID <- mapIds(org.Hs.eg.db,
                             keys=up_GSE150$Gid, #Column containing Ensembl gene ids
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first")
up_GSE150_Entrez <- up_GSE150$EntrezID
up_GSE150_Entrez <- up_GSE150_Entrez[!is.na(up_GSE150_Entrez)]

down_GSE150$EntrezID <- mapIds(org.Hs.eg.db,
                               keys=down_GSE150$Gid, #Column containing Ensembl gene ids
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")
down_GSE150_Entrez <- down_GSE150$EntrezID
down_GSE150_Entrez <- down_GSE150_Entrez[!is.na(down_GSE150_Entrez)]

GO_Up_GSE150 <- basicProfile (up_GSE150_Entrez, onto="BP", level=2, orgPackage="org.Hs.eg.db") 
printProfiles(GO_Up_GSE150, percentage = T)
plotProfiles (GO_Up_GSE150, aTitle="Up-regulated genes GSE150")

GO_Down_GSE150 <- basicProfile (down_GSE150_Entrez, onto="BP", level=2, orgPackage="org.Hs.eg.db") 
printProfiles(GO_Down_GSE150, percentage = T)
plotProfiles (GO_Down_GSE150, percentage = T, aTitle="Down-regulated genes GSE150")

GO_Up_Down_GSE150316 <-mergeProfilesLists(GO_Up_GSE150, GO_Down_GSE150, profNames=c("Up", "Down"))
plotProfiles (GO_Up_Down_GSE150316, percentage=T,aTitle="Up vs Down", legendText = T)

fisherGOProfiles(GO_Up_GSE150$BP, GO_Down_GSE150$BP, method="BH")



