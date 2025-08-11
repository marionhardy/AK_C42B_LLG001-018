
### This script was written by Marion Hardy for a collaboration with 
## Akshaya Karthikeyan from the Lombard lab

## This sets up your DESeq2 object for further data analysis
## This is the C42B cell line treated with Ola at different doses for different 
## lengths

library(tidyverse)
library(readxl)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

## Cut down to the samplkes we want to analyze

data = read.table("./data/rawcountdata_LLG001-027.txt", header = T)
data = data[,1:19]
write.table(data,"./data/rawcountdata_LLG001-018.txt")
rm(data)

## data## Load the raw count matrices, two in this case

counts = read.table("./data/rawcountdata_LLG001-018.txt", header = T)

## Load metadata (sample annotation file)

coldata = read_xlsx("./data/LLG001-018_metadata.xlsx", sheet = "Recoded")
rownames(coldata) = coldata$ID

## Clean up counts matrix and ensembl id

strrep = sub(pattern = "\\.(.*)","",counts$X)
counts$X = strrep
rownames(counts) = counts$X
counts = counts %>% select(!X)

# Create the full model for comparison of samples
# AK said compare DMSO to Ola day 5 and Ola day 9

dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~Treatment) 
dds$Combined = factor(paste0(dds$Set, "_", dds$Treatment))

# checked experimental design, looks good, no recorded confounding variable

# Generate a linear model
dds$Combined = relevel(dds$Combined, "Day1_DMSO")
design(dds) = ~Combined
dds = DESeq(dds)
resultsNames(dds)

# Generate models for the comparisons between mono and bi-therapies
dds$Combined = relevel(dds$Combined, "Day5_DMSO")
design(dds) = ~Combined
dds1 = DESeq(dds)
resultsNames(dds1)

## Checking distribution of counts per sample

as_tibble(assay(dds)) %>%
  gather(sample, value = counts) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample) # looks ok

# Checking size factors and dispersion

sizeFactors(dds) # only takes into account the sequencing depth, looks very clean
plotDispEsts(dds) # verifies normalization, graph looks a-ok

# Checking PCA

rld = vst(dds)

plotPCA(rld,intgroup="Combined") + 
  theme_bw()+
  labs(title = 'PCA per treatment')

# Checking sample similarity

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Combined, sep="-")
colnames(sampleDistMatrix) <- paste(rld$Combined, sep="-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# you can see one of the day9 is clustering with day5
# but the pca shows that very little of the variance is explained
# by the length of Ola treatment so it makes sense that they are similar
# enough to cluster together when computing euclidian distances

## Saving the DESeq object

saveRDS(dds, "./data_output/C42B_DESeq_1-18_Day1_vs_DMSO.Rds")
saveRDS(dds1, "./data_output/C42B_DESeq_1-18_Day5_vs_DMSO.Rds")







