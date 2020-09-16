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
# Use a "base" taxon that is present across all samples at even levels
divnet_asv <- divnet(physeq,
                     base="Firmicutes;Bacilli;Bacillales;Staphylococcaceae;S31;Staphylococcaceae", ncores=4)

### SHANNON AND SIMPSON ALPHA DIVERSITY METRICS AND PLOTS (FIG 1A AND 1B) ###

# Show all estimated alpha and beta diversity metrics
divnet_family %>% names

# Pull them out individually
divnet_family$shannon %>% 
  summary %>%
  add_column("Sample" = physeq %>% otu_table %>% sample_names)

# Extract just the metrics of choice
a <- divnet_family$shannon %>% summary
b <- divnet_family$simpson %>% summary

# Make plots

# First extract just the sample names and metric values
shan <- a %>% select(sample_names, estimate)
simp <- b %>% select(sample_names, estimate)

# Then add Group data (asymptomatic vs symptomatic)
shan <- shan[order(shan$sample_names),]
shan$group <- Noro_baselines_meta$Group

simp <- simp[order(simp$sample_names),]
simp$group <- Noro_baselines_meta$Group

# Plot individual data points

shan %>%
  ggplot(aes(group, estimate, group=group)) +
  geom_point(color="black", size=3, alpha=0.75) +
  theme_classic() # + geom_boxplot()

simp %>%
  ggplot(aes(group, estimate, group=group)) +
  geom_point(color="black", size=3, alpha=0.75) +
  theme_classic() #+ geom_boxplot()

### BETA DIVERSITY: BRAY-CURTIS NMDS PLOT (FIG 1C) ########

# Pull out Bray-Curtis distances as a square distance matrix
bc <- divnet_asv$'bray-curtis'

# Use metaMDS in vegan to run NMDS analysis on Bray-Curtis distance matrix
nmds_bc <- metaMDS(bc, k=3, autotransform=FALSE, engine="isoMDS")

# extract scores from NMDS
data.scores <- as.data.frame(scores(nmds_bc))

# add in metadata (asymptomatic vs symptomatic)
data.scores$Group <- metadata$Group

# Plot nMDS

nmds = data.frame(data.scores, Group=as.factor(data.scores$Group)) 
#shape=as.factor(data.scores$Individual))

p <- ggplot(nmds, aes(x=NMDS1, y=NMDS2, colour=Group)) # shape=shape 

p + theme_classic() + geom_point(size=3) + labs(x = "MDS1", y = "MDS2") +
  scale_color_manual(values=c("red","forestgreen"))

#### RUN ANOSIM ON BRAY-CURTIS DISSIMILARITY MATRIX ####

bc_anosim <- anosim(bc, grouping=metadata$Group)
bc_anosim


#### NMDS and ANOSIM on Mash distance matrix ####

mash <- column_to_rownames(all_k25_full_mash, 'X1')

# Use metaMDS in vegan to run NMDS analysis on Mash distance matrix
nmds_mash <- metaMDS(mash, k=3, autotransform=FALSE, engine="isoMDS")

# extract scores from NMDS
data.scores <- as.data.frame(scores(nmds_mash))

# add in metadata (asymptomatic vs symptomatic)
data.scores$Group <- metadata$Group

# Plot NMDS

nmds = data.frame(data.scores, Group=as.factor(data.scores$Group)) 
#shape=as.factor(data.scores$Individual))

p <- ggplot(nmds, aes(x=NMDS1, y=NMDS2, colour=Group)) # shape=shape 

p + theme_classic() + geom_point(size=3) + labs(x = "MDS1", y = "MDS2") +
  scale_color_manual(values=c("red","forestgreen"))

#### RUN ANOSIM ON MASH DISSIMILARITY MATRIX ####

bc_anosim <- anosim(mash, grouping=metadata$Group)
bc_anosim
