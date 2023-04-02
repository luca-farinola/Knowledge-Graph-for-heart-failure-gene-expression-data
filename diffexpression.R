
# required packages 

require(rstudioapi)
require(dplyr)
require(tidyr)
require(biomaRt)
require(ggplot2)
require(pheatmap)
require(biomaRt)
require(pcaMethods)
require(limma)
require(gridExtra)


# set the directory you are working in 

setwd(dirname(getSourceEditorContext()$path))


# create the directories where i will store important files 
# (file.path() and getwd create reproducible code that can be run on other operating systems)

Imm <- 'immages'
Files <- 'files'

dir.create(file.path(getwd(), Imm ))
dir.create(file.path(getwd(), Files ))

# importing data 

gxData <- read.delim(file = 'data/MAGNET_GeneExpressionData_CPM_19112020.txt', header=TRUE ,sep = '\t',row.names = 1)   
sampleInfo <- read.csv(file = 'data/MAGNET_SampleData_18112022.csv', header=TRUE, sep = ',', stringsAsFactors = TRUE)
geneTotExonLengths <- read.delim("data/MAGNET_exonLengths.txt", as.is = T, row.names = 1) 

# I have two na in diabetes since I want to use Diabetes in the linear model i remove it   

doutlier <- which(is.na(sampleInfo$Diabetes))
sampleInfo <- sampleInfo[-doutlier,]
gxData <- gxData[,-doutlier]

#-----------------------------------------------------------------------------#
# Data exploratoration (genes)                                                #        
#-----------------------------------------------------------------------------#

# I decided to use the fpkm to have to better compare genes as fpkm take into consideration also gene lenght 
# calculation of fpkm 

all(rownames(geneTotExonLengths) == rownames(gxData)) 
cpm2fpkm <- function(x) {
  .t <- 2^(x) * 1E3 / geneTotExonLengths[, 1]
}
gxData_fpkm <- cpm2fpkm(gxData)

# check if col names of the sample are in the same order in metadata file
# table(rownames(t(gxData_fpkm )) == sampleInfo$sample_name)

# add etiology info in gxData_fpkm_et 

gxData_fpkm_et <- gxData_fpkm  
colnames(gxData_fpkm_et) <- sampleInfo$etiology

# aggregate etiology columns toghether 

gxData_fpkm_aggregated <- t(rowsum(t(gxData_fpkm_et), names(gxData_fpkm_et))/c(table(names(gxData_fpkm_et))))

# select the first 100 more variant genes and plot with pheatmap 

var <- apply(gxData_fpkm_aggregated  , 1, var)
selectedGenes <- names(var[order(var, decreasing = T)][1:100])

pdf(file.path(getwd(), Imm,'pheatmap.pdf'),width = 12, height = 6) 
pheatmap(gxData_fpkm_aggregated[selectedGenes,], scale = 'row', show_rownames = FALSE)
dev.off()

#performing PCA with pcaMethods

pcaRes <- pca(t(gxData), nPcs = 10)
plotPCA <- cbind(data.frame(pcaRes@scores), sampleInfo)

# Create plots with the first two principal components 

pdf(file.path(getwd(), Imm,'pca.pdf'),width = 12, height = 6)
ggplot(plotPCA, aes(x = PC1, y = PC2)) + 
  geom_point(aes(col = etiology)) + theme_bw()
dev.off()


#-----------------------------------------------------------------------------#
# Differential expression analysis                                            #        
#-----------------------------------------------------------------------------#

# differential expression with limma i'm taking into consideration age gender and diabetes 
# age is relevant as patients that are 60 years or older are all non failure (NF), additionally PPCM patients are all younger 
# gender both for the imbalance of men and women in the different categories and especially for the contrast ppcm vs nf 
# diabetes as diabetes is ofìbserved only in NF and DCM (with a single exception in HCM) this mean that especially hcm vs nf and ppcm vs nf need for correction 

design <- model.matrix(~0 + etiology + gender + age + Diabetes , data = sampleInfo)
colSums(design)

contrast <- makeContrasts(DCMvsControl = etiologyDCM - etiologyNF,
                          HCMvsControl = etiologyHCM - etiologyNF,
                          PPCMvsControl = etiologyPPCM - etiologyNF, 
                          levels = design)

fit <- lmFit(gxData, design)
fit <- contrasts.fit(fit, contrast)

fit <- eBayes(fit)

results_DCM <- decideTests(fit)
summary(results_DCM)



extract_differentially_expressed <- function(my_contrast, FC, pval) {
  
  Res_table <- topTable(fit, coef = my_contrast, number = nrow(gxData))
  
  # establishing witch gene is invariant and witch is UP/DOWN regulated according to the chosen thresholds 
  
  Res_table$diffexpressed <- "INVARIANT"
  Res_table$diffexpressed[Res_table$logFC > FC & Res_table$adj.P.Val < pval] <- "UP"
  Res_table$diffexpressed[Res_table$logFC < -FC & Res_table$adj.P.Val < pval] <- "DOWN"
  
  # taking the log10 of adjusted p values for visualizzations pourposes and creating the volcano plot 
  
  log10padj <- -log10(Res_table$adj.P.Val)
  
  volcano_plot <- ggplot(data = as.data.frame(Res_table), aes(x = log10padj, y = logFC, color = diffexpressed)) + 
    geom_point() +
    scale_color_manual(values = c("DOWN" = "red", "INVARIANT" = "grey", "UP" = "green")) +
    geom_vline(xintercept=-log10(pval), col="black", linetype = "dotted") +
    geom_hline(yintercept=c(-FC, FC), col="black", linetype = "dotted") + 
    ylim(-4,4) +
    xlim(0,30) +
    theme_bw()
  
  # return both the table of differentially expressed and the volcano plot 
  
  return(list(Res_table[Res_table$diffexpressed != 'INVARIANT',],volcano_plot))
  
  # the output list contains [1] = table of differentialy expressed, [2] = volcanoplot,
}

#choose a teshold I have found 1 and 0.05 to filter for a reasonable ammount of genes  

FC <- 1
pval <- 0.05

# create the result tables and plot the volcanos(second element of the list from extract_differentially_expressed) 

dcmvscontroll <- extract_differentially_expressed('DCMvsControl', FC, pval)

pdf(file.path(getwd(), Imm,'volcano_dcmvscontroll.pdf'),width = 12, height = 6)
dcmvscontroll[2] 
dev.off()

#-----------------------------------------------------------------------------#
# Annotation                                                                  #        
#-----------------------------------------------------------------------------#

#the input of this function is the output of the extract_differentially_expressed (the first element of the list)
# this function gives you the result table with gene symbols column 

getsymbols <- function(genetable) {
  
  genelist <- as.data.frame(genetable[1])
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  genes <- rownames(genelist)
  
  sym_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_id"),values=genes,mart= mart)
  
  genelist <- genelist[sym_list$ensembl_gene_id,]
  
  genelist$hgnc_id <- sym_list$hgnc_id
  
  return(genelist)
  
}

dcmvscontroll_genelist <- getsymbols(dcmvscontroll)


#-----------------------------------------------------------------------------#
# Data export                                                                 #        
#-----------------------------------------------------------------------------#

# save as txt in the file folder 

write.table(dcmvscontroll_genelist, file.path(getwd(), 'Files','dcmvscontroll_genelist.txt'), sep = "\t", 
            row.names = T)


