library(ggplot2)

args <- commandArgs(TRUE)

a <- read.csv(args[1], header=TRUE, sep="\t")
b <- a[, c(1,2)]

colors <- ifelse(b$args[2]>=3.0, "Up", ifelse(b$args[2]<=-3.0, "Down", "Normal"))
b$Type <- colors

jpeg(file=paste(args[2], ".Z.jpeg", sep=""), width=2000, height=1000)
p <- ggplot(b, aes(x=X, y=args[2], color=type)) + geom_point()
p + theme(axis.text.x=element_text(size=10, vjust=0.5, hjust=0.5, angle=90)) +
	labs(x="Location ID") +
	scale_y_continuous(breaks=c(-5.0, -4.5, -4.0, -3.5, -3.0, 3.0, 3.5, 4.0, 4.5, 5.0)) +
	geom_rect(mapping=aes(ymin=-3.0, ymax=3.0, xmin=0, xmax=Inf), alpha=0.01, color=NA)
dev.off()
