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

### GENERATE PHYLOSEQ OBJECT ######

# import ASV table and ID-to-taxonomy table

microbe <- column_to_rownames(Noro_baselines_16S_all, 'Taxonomy')
tax_table <- column_to_rownames(Noro_baselines_taxonomy, 'Taxonomy')
metadata <- column_to_rownames(Noro_baselines_meta, 'Sample')

# convert all to matrices
microbe <- as.matrix(microbe)
tax_table <- as.matrix(tax_table)
metadata <- as.data.frame(metadata)

# combine into phyloseq object
ASV = otu_table(microbe, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
MET = sample_data(metadata)
physeq = phyloseq(ASV, TAX, MET)

### RUN DIVNET USING PHYLOSEQ OBJECT ######

# Run DivNet at the family level for box / scatter plots
divnet_family <- divnet(tax_glom(physeq, taxrank="Family"), ncores=4) 

# Run DivNet at the ASV level for NMDS
divnet_asv <- divnet(physeq, 
                     base="Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;unidentified", ncores=4)

### SHANNON AND SIMPSON ALPHA DIVERSITY METRICS AND PLOTS (FIG S2A AND S2B) ###

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

# Then add Group data (asymptomatic vs symptomatic)
shan <- shan[order(shan$sample_names),]
shan$group <- Noro_baselines_meta$Group

simp <- simp[order(simp$sample_names),]
simp$group <- Noro_baselines_meta$Group

# Plot as box plots with individual data points

shan %>%
  ggplot(aes(group, estimate, group=group))+ geom_boxplot() +
  geom_point(color="black", size=3, alpha=0.75) +
  theme_classic()

simp %>%
  ggplot(aes(group, estimate, group=group))+ geom_boxplot() +
  geom_point(color="black", size=3, alpha=0.75) +
  theme_classic()

### BETA DIVERSITY: BRAY-CURTIS NMDS PLOT (FIG S2C) ########

# Pull out Bray-Curtis distances as a square distance matrix
bc <- divnet_asv$'bray-curtis'

# Use metaMDS in vegan to run NMDS analysis on Bray-Curtis distance matrix
nmds_bc <- metaMDS(bc)

# extract scores from NMDS
data.scores <- as.data.frame(scores(nmds_bc))

# add in metadata (asymptomatic vs symptomatic)
data.scores$Group <- metadata$Group

# Plot nMDS

nmds = data.frame(data.scores, Group=as.factor(data.scores$Group)) 
#shape=as.factor(data.scores$Individual))

p <- ggplot(nmds, aes(x=NMDS1, y=NMDS2, colour=Group)) # shape=shape 

p + theme_classic() + geom_point(size=3) + labs(x = "MDS1", y = "MDS2") 

#### RUN ANOSIM ON BRAY-CURTIS DISSIMILARITY MATRIX ####

bc_anosim <- anosim(bc, grouping=metadata$Group)
bc_anosim
