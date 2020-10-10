#load
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(vegan)
library(reshape)

# load the mash output as distance matrix
df <- data.frame(baselines_mash_dist, row.names=1)
colnames(df) <- rownames(df)
dist_mat <- as.matrix(df)

# load the metadata
metadata <- data.frame(noro_baselines_metadata)

# calculate NMDS with user-supplied distances, 2 dimensions
nmds <- metaMDS(dist_mat, dist, k=2, trace=FALSE, autotransform=FALSE)
nmds <- metaMDS(dist_mat, dist, k=3, trace=FALSE, autotransform=FALSE) #plot=TRUE
# nmds <- metaMDS(dist_mat, dist, k=2, trace=FALSE) 

MDS1 = nmds$points[,1]
MDS2 = nmds$points[,2]
MDS3 = nmds$points[,3]

NMDS = data.frame(MDS1=MDS1, MDS2=MDS2, Group=factor(metadata$Group), 
                  Ind=factor(metadata$Individual)) # MDS3=MDS3,

g <- ggplot(NMDS, aes(x=MDS1, y=MDS2)) + geom_point(aes(col=Group),size=3) + 
  theme_classic() + 
  labs(title="Baseline Metagenomes: Mash Distances (k=25, MDS3)")

g

ggsave("name")