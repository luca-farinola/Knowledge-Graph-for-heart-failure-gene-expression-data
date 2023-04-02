# retrieving ontologies for dcm but it is enough


library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# set the directory you are working in 

setwd(dirname(getSourceEditorContext()$path))

genes <- read.delim(file = 'files/dcmvscontroll_genelist.txt', header=TRUE ,sep = '\t',row.names = 1)
genes <- rownames(genes)

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
go_list <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_id","go_id", "name_1006", "namespace_1003"),values=genes, mart= mart)

significant_Go_BP <- enrichGO(gene = genes,OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP')

significant_Go_BP <- as.data.frame(significant_Go_BP)

significant_Go_MF <- enrichGO(gene = genes,OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'MF')

significant_Go_MF <- as.data.frame(significant_Go_MF)

significant_Go_CC <- enrichGO(gene = genes,OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'CC')

significant_Go_CC <- as.data.frame(significant_Go_CC)

binded <- rbind(significant_Go_BP,significant_Go_MF)
enriched_GO <- rbind(binded,significant_Go_CC)

go_list <- go_list %>%
              filter(go_id %in% enriched_GO$ID)

colnames(go_list)[2] <- "ID"

merged_df <- merge(go_list, enriched_GO, by = "ID")

write.csv(merged_df,'files/geneontology_DCM.txt', row.names = FALSE)


