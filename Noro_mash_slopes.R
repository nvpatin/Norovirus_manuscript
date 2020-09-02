library(ggplot2)
library(readr)
library(reshape2)

d <- data.frame(noro_mash_T0_all_slopes)

slopes <- melt(d, id=c("Individual","Group"))

p <- ggplot(slopes, aes(variable, value, colour=Group))
p + geom_point(size=2, aes(Group, value, colour=Group)) + 
  geom_violin(aes(Group, value, colour=Group), alpha=0) + theme_classic()
