#load
library(phyloseq)
library(DivNet)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(data.table)
library(textshape)
library(magrittr) # to use %>% as a pipe
library(vegan)
library(doParallel)
library(tidyverse)

# import ASV table (eg "table_IDs_sterivex") and ID-to-taxonomy table

microbe <- column_to_rownames(Noro_timeseries_16S_all, 'Taxonomy')
tax_table <- column_to_rownames(Noro_timecourse_taxonomy, 'Taxonomy')
metadata <- column_to_rownames(Noro_timecourse_metadata, 'Sample')

# convert all to matrices
microbe <- as.matrix(microbe)
tax_table <- as.matrix(tax_table)
metadata <- as.data.frame(metadata)

# combine into phyloseq object
ASV = otu_table(microbe, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
MET = sample_data(metadata)
physeq = phyloseq(ASV, TAX, MET)

# Run DivNet at the family level
divnet_family <- divnet(tax_glom(physeq, taxrank="Family"), ncores=4) #X="Group"

# Show all estimated alpha and beta diversity metrics
divnet_family %>% names

# Pull them out individually
divnet_family$shannon %>% 
  summary %>%
  add_column("Sample" = physeq %>% otu_table %>% sample_names)

# Extract just the metrics of choice (don't get fancy)
a <- divnet_family$shannon %>% summary
b <- divnet_family$simpson %>% summary

# Make box plots

# First extract just the sample names and metric values
shan <- a %>% select(sample_names, estimate)
simp <- b %>% select(sample_names, estimate)

# Then add Group data (time point)
shan <- shan[order(shan$sample_names),]
shan$group <- Noro_timecourse_metadata$Time
shan$individual <- as.factor(Noro_timecourse_metadata$Individual)

simp <- simp[order(simp$sample_names),]
simp$group <- Noro_timecourse_metadata$Time
simp$individual <- as.factor(Noro_timecourse_metadata$Individual)

# Plot as violin plots with individual data points overlaid

shan %>%
  ggplot(aes(group, estimate, group=group, color=individual)) +
  geom_point(size=2) + geom_violin(alpha=0) +
  theme_classic()

simp %>%
  ggplot(aes(group, estimate, group=group, color=individual)) + 
  geom_point(size=2) + geom_violin(alpha=0) +
  theme_classic()
