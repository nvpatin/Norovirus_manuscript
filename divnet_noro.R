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

# import ASV table (eg "table_IDs_sterivex") and ID-to-taxonomy table

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

# Run DivNet
divnet_phylum <- divnet(tax_glom(physeq, taxrank="Phylum"), X="Location", ncores=6)
divnet_genus <- divnet(tax_glom(physeq, taxrank="Genus"), X="Group", ncores=4)
divnet_family <- divnet(tax_glom(physeq, taxrank="Family"), X="Group", ncores=4)
divnet_asv <- divnet(physeq, base="580c721cd115033cebadc0d909e4cd83", ncores=4)
divnet_asv <- divnet(physeq, base="56f3f29aaa149dc5b511e5dd1b4c184c", ncores=4)
divnet_asv <- divnet(physeq, base="29cde7953b89ed299cf71fc3632c9745", ncores=4)

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

# Then add time point data
shan <- shan[order(shan$sample_names),]
shan$group <- Noro_baselines_meta$Group

simp <- simp[order(simp$sample_names),]
simp$group <- Noro_baselines_meta$Group

p <- ggplot(shan, aes(group, estimate, group=group))
q <- ggplot(simp, aes(group, estimate, group=group))

q + geom_boxplot() + theme_classic()

# Or plot them
plot(divnet_family$simpson, physeq, col="Group")

divnet_family_group <- physeq %>%
  divnet(X="Group", ncores=4)

plot(divnet_family_time$shannon, physeq, col="Time")

testDiversity(divnet_family, "shannon")

# Pull out Bray-Curtis distances as a square distance matrix
bc <- divnet_family$'bray-curtis'
euc <- divnet_family$'euclidean'

# fix weighted unifrac distance matrix
unw_unifrac <- as.data.frame(distance_matrix)
rownames(unw_unifrac) <- unw_unifrac[,1]
unw_unifrac[,1] <- NULL

PCA <- prcomp(x=bc)
PCAi <- data.frame(PCA$x)

PCA <- prcomp(x=unw_unifrac)
PCAi <- data.frame(PCA$x)

# basic PCA
autoplot(prcomp(x=bc), size=3, 
         label.size=6, frame=FALSE) #data=bc_meta, colour="Individual", 

p <- ggplot(PCAi, aes(x=PC1, y=PC2)) + geom_point(size=3)

# copy data frame and add any metadata columns
bc_meta <- bc
bc_meta$Individual = metadata$Individual
bc_meta$Time = metadata$Time

unw_unifrac_meta <- unw_unifrac
unw_unifrac_meta$Depth = metadata$Depth
unw_unifrac_meta$Season = metadata$Month

PCAi <- data.frame(PCA$x, group=as.factor(bc_meta$Time), 
                   shape=as.factor(bc_meta$Individual))

p <- ggplot(PCAi, aes(x=PC1, y=PC2, col=group, shape=shape)) + geom_point(size=3)

p + theme_classic() + labs(x = "PC1 (87.32%)", y = "PC2 (10.62%)") 

# Use metaMDS in vegan to run NMDS analysis on Bray-Curtis distance matrix
nmds_bc <- metaMDS(bc)
nmds_euc <- metaMDS(euc)

# extract scores from NMDS
data.scores <- as.data.frame(scores(nmds_bc))
data.scores <- as.data.frame(scores(nmds_euc))

# add in metadata
data.scores$Group <- metadata$Group
data.scores$Individual <- Noro_timecourse_metadata$Individual
data.scores$Time <- Noro_timecourse_metadata$Time
data.scores$Month <- BlueHole_all_water_metadata_phyloseq$Month

# Plot nMDS'

nmds = data.frame(data.scores, Group=as.factor(data.scores$Group)) 
           #shape=as.factor(data.scores$Individual))

p <- ggplot(nmds, aes(x=NMDS1, y=NMDS2, colour=Group)) # shape=shape 

p + theme_classic() + geom_point(size=3) + labs(x = "MDS1", y = "MDS2") 


