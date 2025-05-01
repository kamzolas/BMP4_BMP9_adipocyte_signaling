
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("viper")
library(viper)
#install.packages("purrr")
library("purrr")
#install.packages("dplyr")
library("dplyr")



df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}



my_TFA <- function(mypath)
{
  Regulon_file<- read.csv("~/Desktop/TFAs/mouse_network.csv", header=T) ###Open the mouse regulon file
    
  ###subset to the threshold - keep only the most confident TFs
  Regulon_file<- Regulon_file[Regulon_file$confidence=='A'| Regulon_file$confidence=='B'| Regulon_file$confidence=='C' | Regulon_file$confidence=='D',]
  
  
  file.names <- dir(mypath, pattern = "_Degs_with_gene_names.txt")
  for (i in 1:length(file.names)) #Repeat for all the DEGs files
  {
    DEsignature <- read.csv(file = paste0(mypath,file.names[i]), sep = "\t", header = TRUE) #Read the DEGs file
    # Exclude probes with unknown or duplicated gene external_gene_name
    DEsignature = subset(DEsignature, external_gene_name != "" )
    DEsignature = subset(DEsignature, ! duplicated(external_gene_name))
    # Estimatez-score values for the GES. Check VIPER manual for details
    myStatistics = matrix(DEsignature$log2FoldChange, dimnames = list(DEsignature$external_gene_name, 'log2FoldChange') )
    myPvalue = matrix(DEsignature$padj, dimnames = list(DEsignature$external_gene_name, 'padj') )
    mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
    mySignature = mySignature[order(mySignature, decreasing = T)]
    # Estimate TF activities
    mrs = msviper(ges = mySignature, regulon = df2regulon(Regulon_file), ges.filter = F, minsize = 4)
    
    
    TF_activities = data.frame(Regulon = names(mrs$es$nes),
                               Size = mrs$es$size[ names(mrs$es$nes) ], 
                               NES = mrs$es$nes, 
                               p.value = mrs$es$p.value, 
                               FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
    TF_activities = TF_activities[ order(TF_activities$p.value), ]
    # Save results
    write.csv(TF_activities, file = paste0(file.names[i], "_TFAs.csv"))
  }
}



#############################
#############################
#Viper - Transcription Factor Activities using the DEGs obtained by DESeq2
#############################
#############################


#pbat
mypath = "../DEGs pbat/" #Absolute Path with the .txt files (files with the differential expression analysis)
my_TFA(mypath)

#3t3-l1
mypath = "../DEGs 3t3l1/" #Absolute Path with the .txt files (files with the differential expression analysis)
my_TFA(mypath)

  

