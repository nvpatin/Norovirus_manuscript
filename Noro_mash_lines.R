library(ggplot2)

df <- data.frame(noro_mash_T0)
line_types <- c("Mash Distance"=1, "NoV Titer"=2)

p <- ggplot(data=df, aes(x=Day, y=1000*Mash_dist, group=Individual)) + 
  geom_line(aes(color=factor(Individual), linetype="Mash Distance")) + 
  theme_classic() + 
  scale_x_continuous(breaks=c(0,3,6,9,12,15,18,21,24,27,30,33), expand = c(0, 0))

p <- p + geom_line(data=df, aes(x=Day, y=log(NoV_titer), 
                                color=factor(Individual), linetype="NoV Titer"))
                                                                       
p <- p + scale_y_continuous(sec.axis = sec_axis(~exp(.), name="NoV Titer", 
                                                breaks=c(0,1e4,1e6,1e9)), 
                            name="Mash Distance from T0", expand = c(0, 0))

p <- p + labs(color = "Individual", linetype="Data Set") 

p <- p + scale_linetype_manual(values=line_types)

p
