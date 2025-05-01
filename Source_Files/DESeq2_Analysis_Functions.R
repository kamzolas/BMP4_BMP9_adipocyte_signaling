#Differential Expression Analysis using the counts obtained by IRAP (mapper: hisat2, quant_method= htseq2)
#D.E.Analysis_Method=DESeq2
#Functions for all the collaborators

###Libraries and packages
#####
#source("https://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
library("preprocessCore")
#biocLite("sva")
library("sva")
library(RColorBrewer)
library(ggplot2)
library(reshape2)
#biocLite("DESeq2")
#BiocManager::install("DESeq2")
library(DESeq2)
#BiocManager::install("biomaRt")
library('biomaRt')
library('limma')
library("ggrepel")
#####


myDeseq2 <- function(counts_, Control_group, Diet_group, deg_filename, species_dataset)
{
  cts = data.matrix(counts_)
  cts = cts[,c(Control_group,Diet_group)]
  cts <- cts[(rowSums(cts)>dim(cts)[2]),] #Exlude low expressed counts
  
  coldata=matrix("", nrow=dim(cts)[2],ncol=2)
  colnames(coldata) = c("Sample_ID","Healthy_or_NAFLD")
  coldata[,1]=colnames(cts)
  coldata[1:length(Control_group),2]="Healthy"
  coldata[(length(Control_group)+1):(length(Control_group)+length(Diet_group)),2]="NAFLD"
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ Healthy_or_NAFLD) 
  
  dds <- DESeq(dds)
  #res <- results(dds) # A positive Log2FC means that the gene is overexpressed in NAFLD. Respectively, a negative value means that it is overexpressed in Healthy.
  res <- results(dds, cooksCutoff=FALSE) #There were some pvals = NA ... That happened because some genes were outliers and didn't fit the negative bimodal distribution. Deactivating the Cook_Filtering we don't have NA p-values 
  write.table(res, sep = "\t", file =paste(deg_filename,".txt",sep = ""))
  
  #myBiomart2(deg_filename, species_dataset)
  myBiomart(deg_filename, species_dataset)
  myDeseq2Info(deg_filename, counts_, Control_group, Diet_group)
}



myDeseq2_withBatchCorrection <- function(counts_, Control_group, Diet_group, batches, deg_filename, species_dataset)
{
  cts = data.matrix(counts_)
  cts = cts[,c(Control_group,Diet_group)]
  cts <- cts[(rowSums(cts)>dim(cts)[2]),] #Exlude low expressed counts
  
  coldata=matrix("", nrow=dim(cts)[2],ncol=3)
  colnames(coldata) = c("Sample_ID","Healthy_or_NAFLD", "batch_id")
  coldata[,1]=colnames(cts)
  coldata[1:length(Control_group),2]="Healthy"
  coldata[(length(Control_group)+1):(length(Control_group)+length(Diet_group)),2]="NAFLD"
  coldata[,3]= batches[c(Control_group, Diet_group)] # batches is a column with an indication of which sample belongs to which batch. It is a vector with the length of the number of samples
  
  if (length(unique(coldata[,3]))>1)
  {
    print("Run DESeq2 after adjusting for batch effects")
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~batch_id + Healthy_or_NAFLD)
    
    dds <- DESeq(dds)
    #res <- results(dds) # A positive Log2FC means that the gene is overexpressed in NAFLD. Respectively, a negative value means that it is overexpressed in Healthy.
    res <- results(dds, cooksCutoff=FALSE) #There were some pvals = NA ... That happened because some genes were outliers and didn't fit the negative bimodal distribution. Deactivating the Cook_Filtering we don't have NA p-values 
    write.table(res, sep = "\t", file =paste(deg_filename,".txt",sep = ""))
    
    #myBiomart2(deg_filename, species_dataset)
    myBiomart(deg_filename, species_dataset)
    myDeseq2Info(deg_filename, counts_, Control_group, Diet_group)
  }
  else
  {
    print("Run DESeq2 without adjusting for batch effects")
    myDeseq2(counts_, Control_group, Diet_group, deg_filename, species_dataset)
  }
}










myBiomart <- function(deg_filename, species_dataset)
{
  res <- read.table(paste(deg_filename,".txt",sep = ""))
  res = res[,-c(1,3,4)]
  res <- res[order(rownames(res)),]
  mart <- useDataset(species_dataset, useMart("ensembl")) #species_dataset= "rnorvegicus_gene_ensembl" or species_dataset= "mmusculus_gene_ensembl" or species_dataset= "hsapiens_gene_ensembl"
  
  genes=rownames(res)
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                            "external_gene_name"),values=genes,mart= mart)
  res[,4]=rownames(res)
  res_with_gene_names = merge(res,G_list,by.x="V4",by.y="ensembl_gene_id")
  write.table(res_with_gene_names,file=paste(deg_filename,"_with_gene_names.txt",sep = ""),row.names=T,quote=F,sep="\t")
}




myBiomart2 <- function(deg_filename, species_dataset)
{
  #BiocManager::install("org.Mm.eg.db")
  #BiocManager::install("org.Rn.eg.db")
  library(org.Hs.eg.db) #Human
  library(org.Mm.eg.db) #Mouse
  library(org.Rn.eg.db) #Rat
  
  res <- read.table(paste(deg_filename,".txt",sep = ""))

  res[,7] = rownames(res)
  if (species_dataset == "hsapiens_gene_ensembl")
    bm <- AnnotationDbi::select(org.Hs.eg.db,
                                key=res$V7, 
                                columns="SYMBOL",
                                keytype="ENSEMBL")
  if (species_dataset == "mmusculus_gene_ensembl")
    bm <- AnnotationDbi::select(org.Mm.eg.db,
                              key=res$V7, 
                              columns="SYMBOL",
                              keytype="ENSEMBL")
  if (species_dataset == "rnorvegicus_gene_ensembl")
    bm <- AnnotationDbi::select(org.Rn.eg.db,
                                key=res$V7, 
                                columns="SYMBOL",
                                keytype="ENSEMBL")
  library(tidyverse)
  bm <- as_tibble(bm)
  res[,7] = rownames(res)
  res = res[,c(7,1:6)]
  res <- inner_join(res, bm, by=c("V7"= "ENSEMBL"))
  
  write.table(res,file=paste(deg_filename,"_with_gene_names.txt",sep = ""),row.names=T,quote=F,sep="\t")
}





myDeseq2Info <- function(deg_filename, counts_, Control_group, Diet_group)
{
  res <- read.table(file=paste(deg_filename,".txt",sep = ""))

  res = res[!is.na(res$padj),]
  
  statistics = paste("\n\nControl Animals: ", length(Control_group), "\n", toString(colnames(counts_[,Control_group])),
                     "\n\nDiet Animals: ", length(Diet_group), "\n", toString(colnames(counts_[,Diet_group])),
                     "\n\nTotal number of genes after excluding low counts: ", dim(res)[1],
                     "\n\nP_value: ", sum(res$pvalue<.05), " (<0.05)",
                     "\nP_value: ", sum(res$pvalue<.01), " (<0.01)",
                     "\nAdj_P_value: ", sum(res$padj<.05), " (<0.05)",
                     "\nAdj_P_value: ", sum(res$padj<.01), " (<0.01)",
                     "\n\nUpregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange>1), "  (logFC>1 - padj<0.05)",
                     "\nUpregulated Genes: ", sum(res$padj<.01 & res$log2FoldChange>1), "  (logFC>1 - padj<0.01)",
                     "\nUpregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange>2), "  (logFC>2 - padj<0.05)",
                     "\nUpregulated Genes: ", sum(res$padj<.01 & res$log2FoldChange>2), "  (logFC>2 - padj<0.01)",
                     "\nUpregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange>3), "  (logFC>3 - padj<0.05)",
                     "\n\nDownregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange< -1), "  (logFC<-1 - padj<0.05)",
                     "\nDownregulated Genes: ", sum(res$padj<.01 & res$log2FoldChange< -1), "  (logFC<-1 - padj<0.01)",
                     "\nDownregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange< -2), "  (logFC<-2 - padj<0.05)",
                     "\nDownregulated Genes: ", sum(res$padj<.01 & res$log2FoldChange< -2), "  (logFC<-2 - padj<0.01)",
                     "\nDownregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange< -3), "  (logFC<-3 - padj<0.05)",
                     sep = "")
  
  write.table(statistics, sep = "\t", file =paste(deg_filename,"_statistics.txt",sep = ""))
}




New_Deseq2Info <- function(deg_filename, stats_folder, counts_, Control_group, Diet_group)
{
  res <- read.table(file=paste(deg_filename,".txt",sep = ""))
  
  res = res[!is.na(res$padj),]
  
  stats = read.table(file=paste("../",stats_folder,"/SpecificStatistics.txt",sep = ""), header = TRUE)
  
  control_samples = colnames(counts_[,Control_group])
  diet_samples = colnames(counts_[,Diet_group])
  
  rownames(stats) = stats[,1]
  
  statistics = paste("\n\nControl Animals: ", length(Control_group), "\nSample ID: ", toString(stats[control_samples,]$Animal_ID),  "\nGC%: ", toString(stats[control_samples,]$GC.),  "\nReads passed quality filter: ", toString(stats[control_samples,]$Reads_passed_quality_filter),  "\nReads mapped: ", toString(stats[control_samples,]$Reads_mapped),
                     "\n\nDiet Animals: ", length(Diet_group), "\nSample ID: ", toString(stats[diet_samples,]$Animal_ID), "\nGC%: ", toString(stats[diet_samples,]$GC.), "\nReads passed quality filter: ", toString(stats[diet_samples,]$Reads_passed_quality_filter), "\nReads mapped: ", toString(stats[diet_samples,]$Reads_mapped),
                     "\n\nTotal number of genes after excluding low counts: ", dim(res)[1],
                     "\n\nP_value: ", sum(res$pvalue<.05), " (<0.05)",
                     "\nP_value: ", sum(res$pvalue<.01), " (<0.01)",
                     "\nAdj_P_value: ", sum(res$padj<.05), " (<0.05)",
                     "\nAdj_P_value: ", sum(res$padj<.01), " (<0.01)",
                     "\n\nUpregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange>1), "  (logFC>1 - padj<0.05)",
                     "\nUpregulated Genes: ", sum(res$padj<.01 & res$log2FoldChange>1), "  (logFC>1 - padj<0.01)",
                     "\nUpregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange>2), "  (logFC>2 - padj<0.05)",
                     "\nUpregulated Genes: ", sum(res$padj<.01 & res$log2FoldChange>2), "  (logFC>2 - padj<0.01)",
                     "\nUpregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange>3), "  (logFC>3 - padj<0.05)",
                     "\n\nDownregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange< -1), "  (logFC<-1 - padj<0.05)",
                     "\nDownregulated Genes: ", sum(res$padj<.01 & res$log2FoldChange< -1), "  (logFC<-1 - padj<0.01)",
                     "\nDownregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange< -2), "  (logFC<-2 - padj<0.05)",
                     "\nDownregulated Genes: ", sum(res$padj<.01 & res$log2FoldChange< -2), "  (logFC<-2 - padj<0.01)",
                     "\nDownregulated Genes: ", sum(res$padj<.05 & res$log2FoldChange< -3), "  (logFC<-3 - padj<0.05)",
                     sep = "")
  
  write.table(statistics, sep = "\t", file =paste(deg_filename,"_statistics.txt",sep = ""))
}





