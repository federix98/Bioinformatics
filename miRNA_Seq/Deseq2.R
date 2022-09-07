# R version 4.1.3 (2022-03-10) 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

BiocManager::install('EnhancedVolcano')

## https://www.reneshbedre.com/blog/deseq2.html 

# load library 
# DESeq2 version 1.34.0
library("DESeq2")
library(EnhancedVolcano)


install.packages("pheatmap")
library(pheatmap)

install.packages("txtimport")
install.packages("tidyverse")

library(tximport)
library(DESeq2)
library(tidyverse)

setwd("/Users/federico/Desktop/Università/COURSE - Bioinformatics/Progetto/Test")

# get count dataset
# count_matrix <- as.matrix(read.table(file = "/Users/federico/Desktop/Università/COURSE - Bioinformatics/Progetto/Test/Samples-BRMature-miRNAcounts.tsv", sep = '\t', header = TRUE))
# view first two rows
# head(count_matrix, 2)

# get count dataset
count_matrix_Mature <- as.matrix(read.csv(file = "Samples-FULL_Mature-miRNAcounts.csv", sep = ";", row.names = "miRNA"))
# view first two rows
head(count_matrix_Mature, 2)

# At this point the count dataset is loaded
coldata <- data.frame(
  sample = c( "X70459.miRNAcount","X70460.miRNAcount","X70461.miRNAcount","X74665.miRNAcount","X74666.miRNAcount","X74667.miRNAcount","X74668.miRNAcount","X74669.miRNAcount","X74670.miRNAcount","X74671.miRNAcount","X74672.miRNAcount","X74673.miRNAcount","X74674.miRNAcount","X74675.miRNAcount","X74676.miRNAcount","X74677.miRNAcount","X74678.miRNAcount","X74679.miRNAcount","X74680.miRNAcount","X74681.miRNAcount","X74682.miRNAcount","X70462.miRNAcount","X74683.miRNAcount","X74684.miRNAcount" ),
  condition = c( "infected", "infected",  "infected", "infected", "infected", "infected","infected", "infected",  "infected", "infected", "infected", "infected","infected", "infected",  "infected", "infected", "infected", "infected","infected", "infected", "infected","control", "control", "control" ), 
  row.names = "sample" )
coldata$condition <- as.factor(coldata$condition)

all(rownames(coldata) %in% colnames(count_matrix_Mature))

all(rownames(coldata) == colnames(count_matrix_Mature))

dds <- DESeqDataSetFromMatrix(countData = count_matrix_Mature, colData = coldata, 
                              design = ~ condition)

min_counts <- 1

dds <- dds[rowSums(counts(dds)) >= min_counts,]

# set control condition as reference
dds$condition <- relevel(dds$condition, ref = "control")

dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds)

res

res[order(res$padj),]  

write.csv(as.data.frame(res[order(res$padj),] ), file=paste("condition_infected_vs_control_dge_mincounts", min_counts, ".csv", sep = ""))

summary(results(dds, alpha=0.05))


normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)

sample_dists <- assay(dds) %>%
  t() %>%
  dist() %>%
  as.matrix() 

head(sample_dists)

mdsData <- data.frame(cmdscale(sample_dists))
mds <- cbind(mdsData, as.data.frame(colData(dds))) # combine with sample data
head(mds)

rownames(mds)

ggplot(mds, aes(X1, X2, shape = condition, colour = condition,)) + 
  geom_point(size = 3) +
  theme_minimal() + geom_text(aes(label = rownames(mds))) 


genes <- order(res$log2FoldChange, decreasing=TRUE)[1:4]

# more expressed
heatds <- assay(dds[c("hsa-miR-7-5p", "hsa-miR-486-3p", "hsa-miR-133a-3p")])

res[c("hsa-miR-7-5p", "hsa-miR-486-3p", "hsa-miR-133a-3p"), ]

write.csv(res[c("hsa-miR-7-5p", "hsa-miR-486-3p", "hsa-miR-133a-3p"), ], file="mature_expressed.csv")

assay(dds)[genes, ]

samples <- read.csv("samples.csv", sep = ";")
head(samples)



annot_col <- samples %>% column_to_rownames('sample') %>% select(condition) %>% as.data.frame()

head(annot_col)

vsd <- varianceStabilizingTransformation(dds)

pheatmap(assay(dds)[genes, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=annot_col)

colnames(dds)

annot_col
pheatmap(heatds, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=TRUE, annotation_col=annot_col)

EnhancedVolcano(res,
                title = "Volcano plot - Mature miRNA",
                subtitle = bquote(italic(EnhancedVolcano)),
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                ylab = bquote(~-Log[10] ~ italic(P-adj)),
                pCutoff = 0.1,
                FCcutoff = 1, # < -1 downregulated, > 1 upregulated
                ylim = c(0,7),
                xlim = c(-6,6)
                )

EnhancedVolcano(res,
                title = "Volcano plot - Mature miRNA",
                subtitle = bquote(italic(EnhancedVolcano)),
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylab = bquote(~-Log[10] ~ italic(Pvalue)),
)




plotMA(res)





## Genome

# get count dataset
count_matrix_Genome <- as.matrix(read.csv(file = "Samples-FULL_Genome-miRNAcounts.csv", sep = ";", row.names = "miRNA"))
# view first two rows
head(count_matrix_Genome, 2)

# At this point the count dataset is loaded
coldata <- data.frame(
  sample = c( "X70459.genome.miRNAcount","X70460.genome.miRNAcount","X70461.genome.miRNAcount","X74665.genome.miRNAcount","X74666.genome.miRNAcount","X74667.genome.miRNAcount","X74668.genome.miRNAcount","X74669.genome.miRNAcount","X74670.genome.miRNAcount","X74671.genome.miRNAcount","X74672.genome.miRNAcount","X74673.genome.miRNAcount","X74674.genome.miRNAcount","X74675.genome.miRNAcount","X74676.genome.miRNAcount","X74677.genome.miRNAcount","X74678.genome.miRNAcount","X74679.genome.miRNAcount","X74680.genome.miRNAcount","X74681.genome.miRNAcount","X74682.genome.miRNAcount","X70462.genome.miRNAcount","X74683.genome.miRNAcount","X74684.genome.miRNAcount" ),
  condition = c( "infected", "infected",  "infected", "infected", "infected", "infected","infected", "infected",  "infected", "infected", "infected", "infected","infected", "infected",  "infected", "infected", "infected", "infected","infected", "infected", "infected","control", "control", "control" ), 
  row.names = "sample" )
coldata$condition <- as.factor(coldata$condition)

all(rownames(coldata) %in% colnames(count_matrix_Genome))

all(rownames(coldata) == colnames(count_matrix_Genome))

dds <- DESeqDataSetFromMatrix(countData = count_matrix_Genome, colData = coldata, 
                              design = ~ condition)

min_counts <- 1

dds <- dds[rowSums(counts(dds)) >= min_counts,]

# set control condition as reference
dds$condition <- relevel(dds$condition, ref = "control")

dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds)

res

res[order(res$padj),]  

write.csv(as.data.frame(res[order(res$padj),] ), file=paste("condition_infected_vs_control_dge_mincounts_genome_", min_counts, ".csv", sep = ""))

summary(results(dds, alpha=0.05))


normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)

samples <- read.csv("samples_gen.csv", sep = ";")
head(samples)

colnames(dds)

genes <- order(res$log2FoldChange, decreasing=TRUE)[1:10]
assay(dds)[genes, ]

annot_col <- samples %>% column_to_rownames('sample') %>% select(condition) %>% as.data.frame()

write.csv(res[order(res$log2FoldChange, decreasing=TRUE)[1:10],], file="genome_expressed.csv")

dds[genes, ]

pheatmap(assay(dds)[genes, ], cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=TRUE, annotation_col=annot_col)

# more expressed
heatds_gen <- assay(dds[c("hsa-miR-7-5p", "hsa-miR-486-3p", "hsa-miR-133a-3p")])

pheatmap(heatds_gen)

sample_dists <- assay(dds) %>%
  t() %>%
  dist() %>%
  as.matrix() 

head(sample_dists)

mdsData <- data.frame(cmdscale(sample_dists))
mds <- cbind(mdsData, as.data.frame(colData(dds))) # combine with sample data
head(mds)

ggplot(mds, aes(X1, X2, shape = condition, colour = condition,)) + 
  geom_point(size = 4) +
  theme_minimal() + geom_text(aes(label = rownames(mds))) 


plotMA(res)

EnhancedVolcano(res,
                title = "Volcano plot - Genome miRNA",
                subtitle = bquote(italic(EnhancedVolcano)),
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                ylab = bquote(~-Log[10] ~ italic(P-adj)),
                pCutoff = 0.1,
                FCcutoff = 1, # < -1 downregulated, > 1 upregulated
                ylim = c(0, 1.5))


EnhancedVolcano(res,
                title = "Volcano plot - Mature miRNA",
                subtitle = bquote(italic(EnhancedVolcano)),
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylab = bquote(~-Log[10] ~ italic(Pvalue)),
)
annot_col
