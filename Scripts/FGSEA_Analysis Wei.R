#Pathways analysis for DESeq2 data

library(tidyverse)

my_FGSEA <- function(mypath, organism, database)
{
  file.names <- dir(mypath, pattern = "_Degs.txt")
  for (i in 1:length(file.names)) #Repeat for all the DEGs files
  {
    res <- read.table(file = paste0(mypath,file.names[i]), sep = "", header = TRUE) #Read the DEGs file
    
    #select database
    if(database == "GO")
      pathway_GO_or_KEGG_or_REACTOME = "~/Desktop/FGSEA_Pathway Enrichment Analysis/Human_GO_bp_no_GO_iea_symbol.gmt"
    if (database == "KEGG")
      pathway_GO_or_KEGG_or_REACTOME = "~/Desktop/FGSEA_Pathway Enrichment Analysis/c2.cp.kegg.v6.2.symbols.gmt"
    if (database == "REACTOME")
      pathway_GO_or_KEGG_or_REACTOME = "~/Desktop/FGSEA_Pathway Enrichment Analysis/c2.cp.reactome.v6.2.symbols.gmt"
    
    
    #select organism (human or mouse or rat)
    joinby = "ensembl_gene_id"
    if(organism == "HUMAN")
    {
      library(org.Hs.eg.db)
      res[,7] = rownames(res)
      bm <- AnnotationDbi::select(org.Hs.eg.db,
                                  key=res$V7, 
                                  columns="SYMBOL",
                                  keytype="ENSEMBL")
      bm <- as_tibble(bm)
      joinby = "ENSEMBL"
    }
    if (organism == "MOUSE")
      bm <- read.table(file = "~/Desktop/FGSEA_Pathway Enrichment Analysis/mouse_to_human.txt", sep = '\t', header = TRUE) #Read the mouse_to_human table (correspondence using the orthologs)
    if (organism == "RAT")
      bm <- read.table(file = "~/Desktop/FGSEA_Pathway Enrichment Analysis/rat_to_human.txt", sep = '\t', header = TRUE) #Read the rat_to_human table (correspondence using the orthologs)
    
    
    
    
    res[,7] = rownames(res)
    res = res[,c(7,1:6)]
    library(tidyverse)
    res <- inner_join(res, bm, by=c("V7"= joinby))
    
    #Further, all we care about later on is the gene symbol and the test statistic. 
    #Get just those, and remove the NAs. Finally, if we have multiple test statistics for the same symbol, 
    #we’ll want to deal with that in some way. Here I’m just averaging them.
    if(organism == "HUMAN")
    {
      res2 <- res %>% 
        dplyr::select(SYMBOL, stat) %>% 
        na.omit() %>% 
        distinct() %>% 
        group_by(SYMBOL) %>% 
        summarize(stat=mean(stat))
    } else {
     res2 <- res %>% 
      dplyr::select(hsapiens_homolog_associated_gene_name, stat) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(hsapiens_homolog_associated_gene_name) %>% 
      summarize(stat=mean(stat))
    }
    
    #########################################
    #FGSEA
    #########################################
    library(fgsea)
    ranks <- deframe(res2) #The fgsea() function requires a list of gene sets to check, and a named vector of gene-level statistics, where the names should be the same as the gene names in the pathways list. First, let’s create our named vector of test statistics.
    
    #The gmtPathways() function will take a GMT file and turn it into a list. Each element in the list is a character vector of genes in the pathway.
    pathways.hallmark = fgsea::gmtPathways(pathway_GO_or_KEGG_or_REACTOME)
    
    #For KEGG, Lea had another version that had more pathways
    if (database == "KEGG") 
      pathways.hallmark <- readRDS("~/Desktop/FGSEA_Pathway Enrichment Analysis/keggpathwayList_FullList.RDS")
    
    #Now, run the fgsea algorithm with 1000 permutations
    fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500, nperm=1000)
    #sum(fgseaRes[, padj < 0.05])
    
    if (organism == "HUMAN")
    {
      #write.table(as.matrix(fgseaRes)[,1:7],file=paste0("Human_UCAM FGSEA Results/", file.names[i], "_", database,"_PathwayAnalysis.txt"),row.names=T,quote=F,sep="\t")
      #I don't keep the gene names of the pathway. The reason is because I had error. Try including the 8th column and fix the bug!
      # Save as RDS the full fgseaRes matrix
      saveRDS(fgseaRes, file=paste0(file.names[i], "_", database,"_PathwayAnalysis.RDS"))
    }else {
      write.table(as.matrix(fgseaRes)[,1:7],file=paste0(file.names[i], "_", database,"_PathwayAnalysis.txt"),row.names=T,quote=F,sep="\t")
      #I don't keep the gene names of the pathway. The reason is because I had error. Try including the 8th column and fix the bug!
      # Save as RDS the full fgseaRes matrix
      saveRDS(fgseaRes, file=paste0(file.names[i], "_", database,"_PathwayAnalysis.RDS"))
    }
  }
  
}






#############################
#############################
#FGSEA - Pathway Enrichment Analysis using the DEGs obtained by DESeq2
#############################
#############################


  #pbat (mouse)
  mypath = "/Users/kamzolas_macbookpro/Desktop/p-BAT_3T3-L1 (Wei)/Wei_DESeq2/DEGs pbat/"
  for (database in c("KEGG", "REACTOME", "GO"))
    my_FGSEA(mypath, "MOUSE", database)
  
  #3t3-l1 (mouse)
  mypath = "/Users/kamzolas_macbookpro/Desktop/p-BAT_3T3-L1 (Wei)/Wei_DESeq2/DEGs 3t3l1/"
  for (database in c("KEGG", "REACTOME", "GO"))
    my_FGSEA(mypath, "MOUSE", database)

  
  

