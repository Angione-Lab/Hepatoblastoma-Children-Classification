library(matrixStats)
library(tidyr)
library(ggplot2)
library(dplyr)
library(purrr)
library(ggpubr)
library(reshape2)
library(biomaRt)

rm(list=ls()) # clear the environment

##################### change these
data = "FVA_weights" # "genes", "FVA, "pathways" + "_weights", "_clinical_weights", "_all_weights", or "genes_FVA_weights", "FVA_genes_weights", "pathways_genes_weights"
metric = "median_" # "", "median_" 
file_input = paste0("path\\", metric, data, ".csv") # 
file_output = paste0(metric, data, "_barplot.pdf")
#####################

if (metric == "median_") {  ######################## CAMBIARE GITHUB (TUTTO IL BLOCCO IF)
  y_title <- expression("Median weight")
} else {
y_title <- expression("Total weight")
}
dataframe <- read.csv(file = file_input, header = TRUE)   

# convert HGNC IDs into GENE NAME IDs
if (data == "genes_weights" | data == "genes_clinical_weights" | data == "genes_FVA_weights" | data == "genes_all_weights") {
genes <- dataframe[,1]
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
mapping = getBM(attributes=c('hgnc_id', 'hgnc_symbol'), filters = 'hgnc_id', 
                values = genes, mart = ensembl)
names(mapping)[1] <- "Gene"
dataframe <- merge(mapping, dataframe, by="Gene")
}

# to replace the 'unassigned' pathway in the dataframe (pathway data) and set the percentile threshold
if (data == "pathways_weights" | data == "pathways_clinical_weights" | data == "pathways_genes_weights" | data == "pathways_all_weights") {
  dataframe[1,]$Pathway = "unassigned"   # due to alphabetical order it is always in the first row
  percent_threshold = 0.90
} else {
  percent_threshold = 0.995
}

# rename columns
if (data == "FVA_weights" | data == "FVA_clinical_weights" | data == "FVA_genes_weights" | data == "FVA_all_weights") {
  colnames(dataframe) <- c("Omic", "Reaction name", "Weight", "Omic combination")
} else if (data == "pathways_weights" | data == "pathways_clinical_weights" | data == "pathways_genes_weights" | data == "pathways_all_weights") 
{
  colnames(dataframe) <- c("Omic", "Weight", "Omic combination")
} else {
  colnames(dataframe) <- c("Old ID", "Omic", "Weight", "Omic combination")
}

if (data == "pathways_weights" | data == "pathways_clinical_weights" | data == "pathways_genes_weights" | data == "pathways_all_weights") {
  w_save = 5
  h_save = 5
} else if (data == "FVA_weights" | data == "FVA_clinical_weights" | data == "FVA_genes_weights" | data == "FVA_all_weights") {
  w_save = 8.5 
  h_save = 7.5
} else {
  w_save = 5
  h_save = 7.5  
}

p <- dataframe %>% filter(quantile(Weight, percent_threshold)<Weight) %>%
  ggplot(aes(x = reorder(Omic, -Weight), y = Weight)) +
  geom_bar(aes(fill = Weight), stat = "identity", alpha = 0.6) +  
  theme_bw() + 
  xlab("") + ylab(y_title) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        plot.margin=unit(c(0,0,-5,0),"mm"), # to remove white space below
        axis.ticks = element_blank(),
        text = element_text(size=20),   # title on y-axis
        axis.text.x = element_text(hjust = 1, vjust = 1.0, angle = 90, size = 20), # x-axis
        axis.text.y = element_text(size = 20),  # ticks on y-axis
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
  scale_fill_gradient(low = "#E16670", high = "#CC2936") # two colour-codes: (#AF47FF, #6200AC ) for genes, (#E16670, #CC2936) for fluxes

ggsave(filename = file_output, width = w_save, height = h_save, dpi=600, device="pdf")

p
dev.off()
p