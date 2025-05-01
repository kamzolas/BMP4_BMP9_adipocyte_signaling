# Pathway Enrichment Analysis for the mydatagenes using EnrichR


my_EnrichR <- function(geneset, filename_)
{
  my_enriched_paths <- enrichr(geneset, c("KEGG_2019_Mouse", "GO_Molecular_Function_2023", "GO_Biological_Process_2023", "GO_Cellular_Component_2023"))
  write.csv(my_enriched_paths$KEGG_2019_Mouse, file = paste0("KEGG_", filename_, ".csv"))
  write.csv(my_enriched_paths$GO_Molecular_Function_2023, file = paste0("GO_Molecular_Function_", filename_, ".csv"))
  write.csv(my_enriched_paths$GO_Biological_Process_2023, file = paste0("GO_Biological_Process_", filename_, ".csv"))
  write.csv(my_enriched_paths$GO_Cellular_Component_2023, file = paste0("GO_Cellular_Component_", filename_, ".csv"))
  return(my_enriched_paths)
}


library(enrichR)
###find the list of all available databases from Enrichr
dbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Mouse") #Choose the prefered database; GO_Biological_Process_2023 / GO_Cellular_Component_2023 / GO_Molecular_Function_2023 / KEGG_2019_Mouse

mydata = read.csv("../my_data.csv")

pbat = mydata[ ,c(1:3)]
wat = mydata[ ,c(4:6)]

#PBAT
#1. Common in all
#in_all = intersect( intersect(pbat$PBAT_BMP4, pbat$PBAT_BMP9), pbat$PBAT_BMP10)
in_all = intersect(pbat$PBAT_BMP4, pbat$PBAT_BMP9)
in_all = in_all[!is.na(in_all)]
temp = my_EnrichR(in_all, "pbat_common in all BMPs")
#2. only in BMP4
#only_BMP4 = setdiff(pbat$PBAT_BMP4, union(pbat$PBAT_BMP9, pbat$PBAT_BMP10))
only_BMP4 = setdiff(pbat$PBAT_BMP4, pbat$PBAT_BMP9)
only_BMP4 = only_BMP4[!is.na(only_BMP4)]
temp = my_EnrichR(only_BMP4, "pbat_only BMP4")
#3. only in BMP9
#only_BMP9 = setdiff(pbat$PBAT_BMP9, union(pbat$PBAT_BMP4, pbat$PBAT_BMP10))
only_BMP9 = setdiff(pbat$PBAT_BMP9, pbat$PBAT_BMP4)
only_BMP9 = only_BMP9[!is.na(only_BMP9)]
temp = my_EnrichR(only_BMP9, "pbat_only BMP9")

#WAT
#1. Common in all
# in_all = intersect( intersect(wat$T3L1_BMP4, wat$T3L1_BMP9), wat$T3L1_BMP10)
in_all = intersect(wat$T3L1_BMP4, wat$T3L1_BMP9)
in_all = in_all[!is.na(in_all)]
temp = my_EnrichR(in_all, "wat_common in all BMPs")
#2. only in BMP4
#only_BMP4 = setdiff(wat$T3L1_BMP4, union(wat$T3L1_BMP9, wat$T3L1_BMP10))
only_BMP4 = setdiff(wat$T3L1_BMP4, wat$T3L1_BMP9)
only_BMP4 = only_BMP4[!is.na(only_BMP4)]
temp = my_EnrichR(only_BMP4, "wat_only BMP4")
#3. only in BMP9
#only_BMP9 = setdiff(wat$T3L1_BMP9, union(wat$T3L1_BMP4, wat$T3L1_BMP10))
only_BMP9 = setdiff(wat$T3L1_BMP9, wat$T3L1_BMP4)
only_BMP9 = only_BMP9[!is.na(only_BMP9)]
temp = my_EnrichR(only_BMP9, "wat_only BMP9")


#BMP4
#1. Common in both
in_both = intersect(pbat$PBAT_BMP4, wat$T3L1_BMP4)
in_both = in_both[!is.na(in_both)]
temp = my_EnrichR(in_both, "BMP4_in both pbat and wat")
#2. Only in pbat
only_pbat = setdiff(pbat$PBAT_BMP4, wat$T3L1_BMP4)
only_pbat = only_pbat[!is.na(only_pbat)]
temp = my_EnrichR(only_pbat, "BMP4_only pbat")
#3. Only in wat
only_wat = setdiff(wat$T3L1_BMP4, pbat$PBAT_BMP4)
only_wat = only_wat[!is.na(only_wat)]
temp = my_EnrichR(only_wat, "BMP4_only wat")


#BMP9
#1. Common in both
in_both = intersect(pbat$PBAT_BMP9, wat$T3L1_BMP9)
in_both = in_both[!is.na(in_both)]
temp = my_EnrichR(in_both, "BMP9_in both pbat and wat")
#2. Only in pbat
only_pbat = setdiff(pbat$PBAT_BMP9, wat$T3L1_BMP9)
only_pbat = only_pbat[!is.na(only_pbat)]
temp = my_EnrichR(only_pbat, "BMP9_only pbat")
#3. Only in wat
only_wat = setdiff(wat$T3L1_BMP9, pbat$PBAT_BMP9)
only_wat = only_wat[!is.na(only_wat)]
temp = my_EnrichR(only_wat, "BMP9_only wat")

