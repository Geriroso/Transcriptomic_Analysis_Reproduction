###########################################################################
## BALF: R code
###########################################################################

# This is the script containing the whole code for the analysis of the GSE150316 
# dataset.

# Set up paths
wd <- "~/Desktop/TFG/Data/files/CRA002390 and NCBI"
setwd(wd)
data <- file.path(wd, "data")
results <- file.path(wd, "results")
dir.create(data, showWarnings = FALSE)
dir.create(results, showWarnings = FALSE)

# Load packages
library(DESeq2)
library(xlsx)
library(biomaRt)
library(readr)
library(goProfiles)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GenomicFeatures)

# Access the data
f_cnt <- "balf.count.tsv"
cnt_df <- read_delim(file = file.path(data, f_cnt),"\t", escape_double = FALSE, trim_ws = TRUE)
View(cnt_df)

# Treating the raw data
## As I'm using a code extracted from an official study and its processing is
## properly done, I'm going to process the raw data in order to have the genes
## named in 2 different ways and its first processing will be properly done
## Here I'm gonna use the datasets they used to name the genes by their symbol
f_info <- "hg38_gencode.v32.info.tsv"
gene_df <- read_delim(file = file.path(data, f_info), "\t", escape_double = FALSE, trim_ws = TRUE)
gene_info <- gene_df[, c("Gid", "GeneName")]
gene_info <- gene_info[!duplicated(gene_info), ]
names(gene_info) <- c("Gid", "Name")

###########################################################################
## NetworkAnalyst analysis
###########################################################################
df <- data.frame(
  Gid = cnt_df$Gid,
  nCov1 = rowSums(cnt_df[, c("patient1.rep1", "patient1.rep2")]),
  nCov2 = rowSums(cnt_df[, c("patient2.rep1", "patient2.rep2")]),
  Ctrl1 = cnt_df$Ctrl.SRR10571724,
  Ctrl2 = cnt_df$Ctrl.SRR10571730,
  Ctrl3 = cnt_df$Ctrl.SRR10571732)
NAME <- colnames(df)
df <- as.matrix(df)
CLASS_CvsI <- c(rep("Infected", 2), rep("Control", 3))
colnames(df) <- NULL
NetA_BALF <- rbind(NAME, CLASS_CvsI, df)

# As I can't see an easy way to build the needed matrix, I'm gonna add the "#" and the ":" on the .txt file
write.table(NetA_BALF, file = file.path(results, "NetA_BALF_averaged.txt"), sep="\t", quote=F, col.names=NA)
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
f_cnt <- "balf.count.tsv"
cnt_df <- read_delim(file = file.path(data, f_cnt),"\t", escape_double = FALSE, trim_ws = TRUE)

df <- data.frame(
  Gid = cnt_df$Gid,
  nCov1 = rowSums(cnt_df[, c("patient1.rep1", "patient1.rep2")]),
  nCov2 = rowSums(cnt_df[, c("patient2.rep1", "patient2.rep2")]),
  Ctrl1 = cnt_df$Ctrl.SRR10571724,
  Ctrl2 = cnt_df$Ctrl.SRR10571730,
  Ctrl3 = cnt_df$Ctrl.SRR10571732)
rownames(df) <- cnt_df$Gid

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

condition <- factor(c(rep("nCoV", 2), rep("Ctrl", 3)), levels = c("Ctrl", "nCoV"))

dds <- DESeqDataSetFromMatrix(all_cnt_mat, DataFrame(condition), ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res$Gid <- rownames(res)  
res$Symbol <- gene_info$Name[gene_info$Gid %in% rownames(res)]

res$Tag <- "NC"
res$Tag[res$log2FoldChange > log2(fc_cutoff) & res$padj < padj_cutoff] <- "Up"
res$Tag[res$log2FoldChange < -log2(fc_cutoff) & res$padj < padj_cutoff] <- "Down"
res$Tag <- factor(res$Tag, levels = c("Up", "NC", "Down"))

write.table(res, file = file.path(results, "DEA_BALF_NetAStep.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")

sum(res$Tag == "Up")
sum(res$Tag == "Down")  

# Metascape Analysis
## We need to prepare the DEGs in order to introduce them on the Metascape software
## Here we keep the genes that have been tagged depending on their FC and p-adj
up_BALF <- subset(res, res$Tag=="Up")
down_BALF <- subset(res, res$Tag=="Down")
title <- c("Up", "Down")

## We need to apply the requiered format from Metascape
up_down <- matrix(0, ncol = 2, nrow = 2)
up_down <- data.frame(up_down)
up_down[,1] <- title
colnames(up_down) <- c("Names", "Genes")
up_genes <- paste(up_BALF$Symbol, collapse = ", ")
up_down[1,2] <- up_genes
down_genes <- paste(down_BALF$Symbol, collapse = ", ")
up_down[2,2] <- down_genes
View(up_down)
write.table(up_down, file = file.path(results,"Up&Down-regulated_BALF_NetAStep.txt"),sep = "\t", quote = F, row.names = F, col.names = T)

###########################################################################
## BALF DESeq2 analysis: Default settings from the DESeq2 Package
###########################################################################

# Creation of the DESeq object
if(load.data){
  load(file = file.path(data, "Counts.Rda"))
  load(file = file.path(data, "Top.Rda"))
}else{
  cts <- read.table(file = file.path(data, "balf_raw_count.txt"), header = T, sep = "\t", row.names = 1) 
  condition <- factor(c(rep("nCoV", 4), rep("Ctrl", 3)), levels = c("Ctrl", "nCoV"))
  samples <- factor(c(rep(1:5, c(2,2,1,1,1))))
  coldata <- data.frame(condition, samples)
  rownames(coldata) <- c(
    "patient1.rep1",
    "patient1.rep2",
    "patient2.rep1",
    "patient2.rep2",
    "Ctrl.SRR10571724",
    "Ctrl.SRR10571730",
    "Ctrl.SRR10571732"
  )
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
  
  # Pre-filtering to remove rows with low counts
  keep <- rowSums(counts(dds)) >= 10 
  dds <- dds[keep,]
  
  # Collapsing the technical replicates
  dds$run <- paste0("sample",1:7)
  ddsColl <- collapseReplicates(dds, dds$samples, dds$run, renameCols = F)
  colnames(ddsColl[,1:2]) <- paste0("Patient", 1:2) 
  
  # DEA
  ## Here, I conduct the DEA with the collapsed DESeq object
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
write.xlsx(res, file = file.path(results, "DEA_BALF_DS.xlsx"), sheetName = "DESeq2 data")
write.table(res, file = file.path(results, "DEA_BALF_DS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")

#TPM
# Extract gene length
library(GenomicFeatures)

txdb <- makeTxDbFromGFF(file = file.path(data, "gencode.v32.annotation.gtf"),format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
length(exons.list.per.gene)
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
exonic.gene.sizes <- as.data.frame(exonic.gene.sizes)
colnames(exonic.gene.sizes) <- "Gene_length"
exonic.gene.sizes$Gid <- rownames(exonic.gene.sizes)
View(exonic.gene.sizes)
exonic.gene.sizes <- exonic.gene.sizes[, -2]

#Calculate gene length and TPM
cnt_df <- read_delim(file = file.path(data, f_cnt),"\t", escape_double = FALSE, trim_ws = TRUE)
df <- data.frame(
  Gid = cnt_df$Gid,
  nCov1 = rowSums(cnt_df[, c("patient1.rep1", "patient1.rep2")]),
  nCov2 = rowSums(cnt_df[, c("patient2.rep1", "patient2.rep2")]),
  Ctrl1 = cnt_df$Ctrl.SRR10571724,
  Ctrl2 = cnt_df$Ctrl.SRR10571730,
  Ctrl3 = cnt_df$Ctrl.SRR10571732)
df <- merge(df, gene_info, all = FALSE)
sum(duplicated(df_cnt$Name))

# Know the duplicated rows names
n_occur <- data.frame(table(df_gl$Name))
repeated <- n_occur[n_occur$Freq > 1,]
repeated_names <- as.character(repeated$Var1)
df_gl[df_gl$Name %in% n_occur$Var1[n_occur$Freq > 1],]

df_gl <- merge(df, exonic.gene.sizes, all=FALSE)
View(df_gl)

df_cnt <- data.frame(
  Name = df$Name,
  nCov1 = rowSums(cnt_df[, c("patient1.rep1", "patient1.rep2")]),
  nCov2 = rowSums(cnt_df[, c("patient2.rep1", "patient2.rep2")]),
  Ctrl1 = cnt_df$Ctrl.SRR10571724,
  Ctrl2 = cnt_df$Ctrl.SRR10571730,
  Ctrl3 = cnt_df$Ctrl.SRR10571732)
View(df_cnt)

# I eliminate the duplicated rows since they are very similar
df_cnt <- df_cnt[!duplicated(df_cnt$Name),]
df_gl <- df_gl[!duplicated(df_gl$Name),]
rownames(df_cnt) <- df_cnt$Name
df_cnt <- df_cnt[,-1]


tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## En principi ho és, però seria necessari assegurar-se que l'ordre dels gens és el mateix per no barrejar gene lengths

tpms <- apply(df_cnt, 2, function(x) tpm(x, df_gl$Gene_length))
View(tpms)
write.table(tpms, file = file.path(results, "TPM_counts_BALF.txt"), sep = "\t", quote = F, col.names = NA)

###########################################################################
## CPM from all matrix
###########################################################################
cnt_df <- read_delim(file = file.path(data, f_cnt), "\t", escape_double = FALSE, trim_ws = TRUE)
df <- data.frame(
  Gid = cnt_df$Gid,
  nCov1 = rowSums(cnt_df[, c("patient1.rep1", "patient1.rep2")]),
  nCov2 = rowSums(cnt_df[, c("patient2.rep1", "patient2.rep2")]),
  Ctrl1 = cnt_df$Ctrl.SRR10571724,
  Ctrl2 = cnt_df$Ctrl.SRR10571730,
  Ctrl3 = cnt_df$Ctrl.SRR10571732)

cnt <- df[, 2:ncol(df)]
rownames(cnt) <- df$Gid
cnt_mat <- as.matrix(cnt)
View(norm_expr_df)

res_df <- expr_res[, c("Gid", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Tag")]
norm_mat <- 1e6 * t(t(cnt_mat) / colSums(cnt_mat))
norm_expr_df <- as.data.frame(norm_mat)
norm_expr_df$Gid <- rownames(norm_mat)
norm_expr_df <- right_join(gene_info, norm_expr_df)
norm_expr_df <- norm_expr_df[!duplicated(norm_expr_df$Name),]
norm_expr_df <- as.data.frame(norm_expr_df)
rownames(norm_expr_df) <- norm_expr_df$Name
norm_expr_df <- norm_expr_df[, -1]
norm_expr_df <- norm_expr_df[, -1]
write.table(norm_expr_df, file = file.path(results, "CPM_all_counts_BALF.txt"), sep = "\t", quote = F, col.names = NA)

###########################################################################
## Extra
###########################################################################

###########################################################################
## Excel
###########################################################################
Ensembl_version_id <- df$Gid
EnsemblID <- gsub("\\..*","",Ensembl_version_id)

#columns(org.Hs.eg.db) # returns list of available keytypes
EntrezID = mapIds(org.Hs.eg.db,
                  keys= EnsemblID , #Column containing Ensembl gene ids
                  column="ENTREZID",
                  keytype="ENSEMBL",
                  multiVals="first")

Symbols = gene_info$Name

df_cnt <- df[,-1]
final_mat <- cbind(Symbols, EntrezID, EnsemblID, Ensembl_version_id, df_cnt)
rownames(final_mat) <- Ensembl_version_id
final_mat <-merge(final_mat, expr_res, by = 0, all = FALSE)
final_mat <- final_mat[, -1]
write.table(final_mat, file = file.path(results, "BALF_results_matrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ",")

###########################################################################
## GoProfiles
###########################################################################
up_BALF <- subset(expr_res, expr_res$Tag=="Up")
down_BALF <- subset(expr_res, expr_res$Tag=="Down")


up_BALF$EntrezID <- mapIds(org.Hs.eg.db,
                           keys=up_BALF$Symbol, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
up_BALF_Entrez <- up_BALF$EntrezID
up_BALF_Entrez <- up_BALF_Entrez[!is.na(up_BALF_Entrez)]

down_BALF$EntrezID <- mapIds(org.Hs.eg.db,
                             keys=down_BALF$Symbol, #Column containing Ensembl gene ids
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first")
down_BALF_Entrez <- down_BALF$EntrezID
down_BALF_Entrez <- down_BALF_Entrez[!is.na(down_BALF_Entrez)]

GO_Up_BALF <- basicProfile (up_BALF_Entrez, onto="BP", level=2, orgPackage="org.Hs.eg.db") 
printProfiles(GO_Up_BALF, percentage = T)
plotProfiles (GO_Up_BALF, aTitle="Up-regulated genes BALF")

GO_Down_BALF <- basicProfile (down_BALF_Entrez, onto="BP", level=2, orgPackage="org.Hs.eg.db") 
printProfiles(GO_Down_BALF, percentage = T)
plotProfiles (GO_Down_BALF, percentage = T, aTitle="Down-regulated genes BALF")

GO_Up_Down_BALF <-mergeProfilesLists(GO_Up_BALF, GO_Down_BALF, profNames=c("Up", "Down"))
plotProfiles (GO_Up_Down_BALF, percentage=T,aTitle="Up vs Down", legendText = T)

fisherGOProfiles(GO_Up_BALF$BP, GO_Down_BALF$BP, method="BH")

