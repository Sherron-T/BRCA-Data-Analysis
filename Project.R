##Set Location of the Project Folder
setwd("/Users/sherront/Documents/UT 2020-2021/BCH 339N/Final Project/")

male_names_list <- "TCGA-A1-A0SM
TCGA-A8-A085
TCGA-AC-A62V
TCGA-AO-A1KQ
TCGA-AQ-A54O
TCGA-AR-A1AV
TCGA-BH-A0B4
TCGA-BH-A0DD
TCGA-D8-A1XS
TCGA-E2-A14W
TCGA-EW-A1PD
TCGA-EW-A6SA"

female_names_list <- "TCGA-A1-A0SO
TCGA-A2-A0CX
TCGA-A2-A25A
TCGA-A8-A0A7
TCGA-A8-A09B
TCGA-A8-A09N
TCGA-B6-A0I9
TCGA-BH-A202
TCGA-C8-A1HM
TCGA-C8-A8HQ
TCGA-OL-A5RW
TCGA-PE-A5DC"

## Initialize list of male and female pateint IDs
male_names <- strsplit(male_names_list, "\n")[[1]]
female_names <- strsplit(female_names_list, "\n")[[1]]


## Reads the DNA expression file and the read counts from the different patient data files
male_data_vector <- c()
female_data_vector <- c()
for(i in 1:length(male_names)){
  string_name <- paste(c("MALE",male_names[i],male_names[i]),collapse="/")
  file <- read.csv(string_name, sep = "\t", header = FALSE)
  male_data_vector <- c(male_data_vector, file[[2]])
}
for(i in 1:length(female_names)){
  string_name <- paste(c("FEMALE",female_names[i],female_names[i]),collapse="/")
  file <- read.csv(string_name, sep = "\t", header = FALSE)
  female_data_vector <- c(female_data_vector, file[[2]])
}
str(male_data_vector)
str(female_data_vector)
genes_list <- file[[1]]

library(tidyverse)
library(biomaRt)

## Retrieves the gene names from bioMart for the ensembl/gencode IDs
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
simplified_genes_list <- sub("[.][0-9]*","",genes_list) ## Formats from gencode to ensembl
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = sub("[.][0-9]*","",genes_list), mart= mart)
renamed_list <- c()
for(i in 1:length(simplified_genes_list)){
  named_gene <-gene_IDs[gene_IDs$ensembl_gene_id == simplified_genes_list[i],]$hgnc_symbol
  if(length(named_gene)>1){
    named_gene <- paste(named_gene, collapse = "+")
  }
  if(rlang::is_empty(named_gene)){
    renamed_list <- c(renamed_list, simplified_genes_list[i] )
  }
  else{
    renamed_list <- c(renamed_list,named_gene)
  }
}

## Reformats patients IDs to have the M/F abbreviation
## Makes the RNA sequence counts matrix for the DESeq2 Analysis
male_names <- paste("M",male_names,sep="-")
female_names <- paste("F",female_names,sep="-")
rnaseq_counts <- matrix(c(male_data_vector,female_data_vector), nrow=60488, ncol = 24)
colnames(rnaseq_counts) = c(male_names,female_names)
row.names(rnaseq_counts) = renamed_list


library('pheatmap')
## Heatmap of the correlation between the samples for all genes 
cors = cor(rnaseq_counts, method="spearman")
pheatmap(cors, cluster_cols = F, cluster_rows = F)

## Heatmap of the correlation between the samples for the interested genes 
cors = cor(rnaseq_counts[c("BRCA1","BRCA2","PALB2","CHEK2","CDH1","PTEN","STK11","TP53","AR","ESR1","TRH"),], method="spearman")
pheatmap(cors, cluster_cols = F, cluster_rows = F)


## Assigns the factors of sex for the DNA expression analysis
all_names <- c(male_names,female_names)
Factor_level_1 = male_names
Factor_level_2 = female_names
factor_of_interest = all_names %in% Factor_level_1
factor_of_interest = as.factor(factor_of_interest)
factor_of_interest = relevel(factor_of_interest, ref = "TRUE")

library('DESeq2')
library('R.utils')
library('data.table')

## Initializes the DESeq Matrix
dds <- DESeqDataSetFromMatrix (countData = rnaseq_counts,
                               colData = DataFrame(Variable = factor_of_interest ),
                               design = ~ Variable)

counts_per_million <- function (count_matrix) {
  for(column in 1:ncol(count_matrix)){
    library_size <- sum(count_matrix[,column])
    library_size_per_mil <- library_size / 1000000
    count_matrix[,column] <- count_matrix[,column] / library_size_per_mil
  }
  return(count_matrix)
}

## Normalizes the counts for the matrix
rpm <- counts_per_million(rnaseq_counts)

## Only remove the genes with the low number of counts from more than 3 samples
keep <- attr(which(apply(rpm > 1, 1, function(x) sum(x) > 3)),"names")
dds <- dds[keep,]
dds <- DESeq(dds)

## lfcShrink reduces the noise and variability associated from genes with low read counts
resLFC <- lfcShrink(dds, coef= "Variable_FALSE_vs_TRUE")
resLFC <- resLFC[order(resLFC$pvalue),]
summary(resLFC, alpha = 0.05)
results(dds)
plotMA(resLFC)

## Generate a report for the scatterplot graphs of the genes 
library("ReportingTools")
library("apeglm")
des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./report")
publish(dds,des2Report, pvalueCutoff=1,
        annotation.db="org.Hs.eg", factor = factor_of_interest,
        reportDir="./reports")
finish(des2Report)
resLFC

## Subset the genes of interest and make another report for the scatterplots
dds2 <- dds[c("BRCA1","BRCA2","PALB2","CHEK2","CDH1","PTEN","STK11","TP53","AR","ESR1","TRH"),]
publish(dds2,des2Report, pvalueCutoff=1,
        annotation.db="org.Hs.eg", factor = factor_of_interest,
        reportDir="./reports")

short <- resLFC[c("BRCA1","BRCA2","PALB2","CHEK2","CDH1","PTEN","STK11","TP53","AR","ESR1","TRH"),]
write.csv(resLFC[c("BRCA1","BRCA2","PALB2","CHEK2","CDH1","PTEN","STK11","TP53","AR","ESR1","TRH"),],"table.csv")
write.csv(resLFC,"resLFC.csv")

ggplot(data = read.csv("table.csv"),
       aes(x = reorder(X,log2FoldChange), y = log2FoldChange,
           fill = log2FoldChange > 0))+
  geom_bar(stat = "identity")+
  coord_flip()+
  xlab("Gene")+
  theme(legend.position="none")+
  theme(text = element_text(size = 15))  
##negative log2fold = higher male

install.packages("viridis")  # Install
library("viridis") 

ggplot(read.csv("table.csv"), aes(x=reorder(X,-pvalue), fill=as.factor(X), y=pvalue )) + 
  geom_bar( stat='identity') + 
  theme(legend.position="none") +
  scale_fill_viridis(option="magma",discrete = TRUE) +
  xlab("Gene")+
  theme(text = element_text(size = 15))

## Subset genes correlation for the patients
subset_genes = list("BRCA1", "BRCA2", "PALB2", "CHEK2", "CDH1", "PTEN", "STK11",
                    "TP53", "AR", "ER", "TRH")
subset_rnaseq_counts = subset(rnaseq_counts, rownames(rnaseq_counts) %in% subset_genes)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
subset_rnaseq_counts_norm <- t(apply(subset_rnaseq_counts, 1, cal_z_score))
pheatmap(subset_rnaseq_counts_norm, cluster_cols = F)

## Metadata plot
data <- read.csv("Sample Metadata.csv")
theme_set(theme_bw())
ggplot(data, aes(x = case_submitter_id, y=pathologic_stage)) + 
  geom_point(stat='identity', aes(col=gender), size=3)  +
  theme(legend.position="none") + 
  scale_color_manual(name="gender", 
                     labels = c("male", "female"), 
                     values = c("male"="#0064d4", "female"="#ec4cf6")) +
  scale_x_discrete(limits=data$case_submitter_id)+
  ylab("Pathological Stage")+
  xlab("Patient Case ID")+
  coord_flip()