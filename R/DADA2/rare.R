library(ggplot2)

num <- read.table("observed_species.txt", sep="\t", header=T)

ggplot(num,aes(x=OTU_Num, y=OB, colour=Samples)) + geom_line() +
  labs(x="Number of Sequences", y="Observed Species", title="Observed Species Rarefaction Curve", fill="") +
  scale_x_continuous(breaks=c(0, 5000, 10000, 15000, 20000, 25000)) +
  scale_y_continuous(breaks=c(0, 200, 400, 600, 800, 1000, 1200, 1400)) +
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA))+
  scale_fill_discrete(name="Samples") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"))
