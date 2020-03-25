#
library(ggplot2)

data <- read.table("obser_count_1000.txt", sep="\t", header=T)

ggplot(data,aes(x=Num, y=OB, colour=Samples))+
  theme_bw()+
  geom_line()+
  scale_y_continuous(breaks=c(0, 200, 400, 600, 800, 1000, 1200))+
  xlab("Number of Sequences") + ylab("Observed Species")+
  theme(panel.grid.major=element_blank(), panel.grid.minor= element_blank(),
        panel.background=element_blank(), axis.line=element_line(colour="black"))
  

