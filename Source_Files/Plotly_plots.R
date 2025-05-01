
#install.packages("plotly")
library(plotly)

#my_rplotly is a function to plot the interactive Volcano Plots
my_rplotly <- function(deg_filename)
{
  degs <- read.csv(paste(deg_filename,"_with_gene_names.txt",sep = ""), sep = "\t")
  degs_without_names <- read.table(paste(deg_filename,".txt",sep = ""))
  degs_without_names = degs_without_names[ is.element(rownames(degs_without_names),degs$V7),]
  
  log2FoldChange=data.frame(degs$log2FoldChange)
  AdjPvalue=data.frame(degs$padj)
  
  
  
  colors_ = c(rep("adj_pval<0.05",dim(degs)[1]))
  colors_[which(abs(degs$log2FoldChange) >2)] = "|LogFC|>2"
  colors_[which(degs$padj >.05)] = "adj_pval>0.05"
  
  
  p <- plot_ly(x=~log2FoldChange[,], y=~-log10(AdjPvalue)[,], type = 'scatter', mode='markers', color = ~colors_,
        text=paste("Gene:",degs$external_gene_name,"\nLogFC:",formatC(degs$log2FoldChange,digits=2),"\nPvalue:",formatC(degs$pvalue,digits=2),
                   "\nAdj.Pvalue:",formatC(degs$padj,digits=2)))
  p <- p %>% layout(
    title = 'Volcano Plot',
    xaxis = list(
      title = "log2FoldChange",
      range = c( -max(log2FoldChange) -.2, max(log2FoldChange) +.2)
    ),
    yaxis = list(
      title = "-log10AdjPvalue"
    )
  )
  htmlwidgets::saveWidget(as_widget(p), paste(deg_filename,"_Volcano_padj.html")) #Save as html
}





#myPathways is a function to plot the heatmaps with the pathways in the different timepoints
#It prints 2 heatmaps; one for the pathways that have been obtained according to the upregulated genes and one with the pathways that have been obtained according to the down regulated genes

library(pheatmap)
library(enrichR)

myPathways <- function(mypath, gene_padj, gene_log2FC_Up, gene_log2FC_Down, path_padj, dbs, c, Row_size_Up, Row_size_Down)
{
  file.names <- dir(mypath, pattern = "with_gene_names.txt")
  file.names = file.names[c]
  up_pathways = down_pathways = c()
  all_up=all_down=number_of_paths_up=number_of_paths_down=NULL
  info_for_genes_and_paths = c()
  
  for (i in 1:length(file.names)) 
  {
    file_ <- read.table(file = paste(mypath,file.names[i], sep = ""), sep = '\t', header = TRUE)
    up_reg = as.character(file_[which(file_['padj'] < gene_padj & file_['log2FoldChange'] > gene_log2FC_Up),5])  #Thresholds for adj.pvalues and logfoldchange
    down_reg = as.character(file_[which(file_['padj'] < gene_padj & file_['log2FoldChange'] < gene_log2FC_Down),5])
    my_enriched_up <- enrichr(up_reg, dbs)
    my_enriched_down <- enrichr(down_reg, dbs)
    
    paths_up = my_enriched_up$KEGG_2016[,1]
    paths_down = my_enriched_down$KEGG_2016[,1]
    up = my_enriched_up$KEGG_2016 [,c(1,4)]
    down = my_enriched_down$KEGG_2016 [,c(1,4)]
    
    up_significant = up[which(up[,2] < path_padj),] #Significant pathways for Up-regulated genes --> Column1=pathway_names, Column2=Pathway_adj.p_values
    down_significant = down[which(down[,2] < path_padj),]
    all_up = rbind(all_up, up_significant) #Keep all the significant pathways
    number_of_paths_up = rbind(number_of_paths_up, length(up_significant[,1])) #How many paths from each time-point
    all_down = rbind(all_down, down_significant)
    number_of_paths_down = rbind(number_of_paths_down, length(down_significant[,1]))
    up_pathways = union(up_pathways, up_significant[,1])
    down_pathways = union(down_pathways, down_significant[,1])
    
    info_for_genes_and_paths = c(info_for_genes_and_paths, length(up_reg), length(down_reg), length(paths_up), length(paths_down), dim(up_significant)[1],dim(down_significant)[1])
  }
  
  column_names = file.names
  for (j in 1:length(file.names)) 
  {
    column_names[j] = strsplit(file.names,"_")[[j]][1]
  }
  
  
  #Heatmaps for the "Up" Pathways
  m_up <- matrix(0, nrow = length(up_pathways), ncol = length(file.names))
  if(dim(m_up)[1]>2) #We create the heatmaps only if 2 or more pathways have been detected
  {
    rownames(m_up) = up_pathways
    colnames(m_up) = column_names
    
    v=c()
    u=0
    for (i in 1:length(number_of_paths_up)) 
    {
      u=u+1
      v=c(v,(rep(u,number_of_paths_up[i])))
    }
    all_up[,"Group"] = v
    
    for (i in 1:dim(all_up)[1])
    {
      pathway = all_up[i,1]
      for (j in 1:dim(m_up)[1]) 
      {
        if(pathway == rownames(m_up)[j])
        {
          m_up[pathway,all_up[i,3]] = round(-log10(all_up[i,2]),2)
        }
      }
    }
    
    #plot_ly(x = colnames(m_up), y = rownames(m_up),
    #        z = m_up, type = "heatmap", colors = colorRamp(c("white", "red")))
    pdf(file=paste("UpPathways_GeneFC>",gene_log2FC_Up,"_Gene_padj<", gene_padj,"_Pathway_padj<", path_padj,".pdf"))
    #pheatmap(m_up,trace="none", col=colorRampPalette(c("white","red"),bias=1)(101),cluster_rows=T, cluster_cols=F)
    pheatmap(m_up,trace="none", col=colorRampPalette(c("white","red"),bias=1)(101),cluster_rows=T, cluster_cols=F, fontsize_row = Row_size_Up)   
    dev.off()
  }
  
  
  #Heatmaps for the "Down" Pathways
  m_down <- matrix(0, nrow = length(down_pathways), ncol = length(file.names))
  if(dim(m_down)[1]>2) #We create the heatmaps only if 2 or more pathways have been detected
  {
    rownames(m_down) = down_pathways
    colnames(m_down) = column_names
    
    v=c()
    u=0
    for (i in 1:length(number_of_paths_down)) 
    {
      u=u+1
      v=c(v,(rep(u,number_of_paths_down[i])))
    }
    all_down[,"Group"] = v
    
    for (i in 1:dim(all_down)[1])
    {
      pathway = all_down[i,1]
      for (j in 1:dim(m_down)[1]) 
      {
        if(pathway == rownames(m_down)[j])
        {
          m_down[pathway,all_down[i,3]] = round(-log10(all_down[i,2]),2)
        }
      }
    }
    
    #plot_ly(x = colnames(m_down), y = rownames(m_down),
    #        z = m_down, type = "heatmap", colors = colorRamp(c("white", "red")))
    pdf(file=paste("DownPathways_GeneFC<",gene_log2FC_Down,"_Gene_padj<", gene_padj,"_Pathway_padj<", path_padj,".pdf"))
    pheatmap(m_down,trace="none", col=colorRampPalette(c("white","red"),bias=1)(101),cluster_rows=T, cluster_cols=F, fontsize_row = Row_size_Down)
    dev.off()
  }
  myPathwaysInfo(file.names, gene_padj, gene_log2FC_Up, gene_log2FC_Down, path_padj, dbs, info_for_genes_and_paths, m_up, m_down)
}



myPathwaysInfo <- function(file.names, gene_padj, gene_log2FC_Up, gene_log2FC_Down, path_padj, dbs, info_for_genes_and_paths, m_up, m_down)
{
  statistics = paste("\n\n\nStatisticts for the Pathway Enrichment Analysis using EnrichR\nAdj_P_value threashold for the genes (up and down regulated): ",gene_padj,
                     "\nLog2FC for the up-regulated genes: ",gene_log2FC_Up,
                     "\nLog2FC for the down-regulated genes: ",gene_log2FC_Down,
                     "\nAdj_P_value for the significant pathways: ",path_padj,
                     "\nEnrichR database for our analysis: ",dbs, 
                     sep = "")
  
  for (i in 1:length(file.names)) 
  {
    statistics = paste(statistics, "\n\nFile", i, ": ", file.names[i],
                       "\nUp-regulated genes: ", info_for_genes_and_paths[(i-1)*6+1],
                       "\nDown-regulated genes: ",info_for_genes_and_paths[(i-1)*6+2],
                       "\nPathways for Up-regulated genes: ", info_for_genes_and_paths[(i-1)*6+3], "  (significants: ", info_for_genes_and_paths[(i-1)*6+5],")", 
                       "\nPathways for Down-regulated genes: ", info_for_genes_and_paths[(i-1)*6+4],"  (significants: ", info_for_genes_and_paths[(i-1)*6+6],")",
                       sep="")
    
  }
  write.table(statistics, sep = "\t", file = paste("Pathways_GeneFC>",gene_log2FC_Up,"_<",gene_log2FC_Down,"_Gene_padj<",gene_padj,"_Path_padj<",path_padj,"statistics.txt"))
}
