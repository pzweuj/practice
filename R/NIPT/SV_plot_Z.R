library(ggplot2)
args <- commandArgs(TRUE)


inputfile <- read.csv(args[1], header=TRUE, sep="\t")
data_use <- inputfile[, c(1, 2)]
names(data_use) <- c("ID", "Data")

colors <- ifelse(data_use[2]>=3.0, "Up", ifelse(data_use[2]<=-3.0, "Down", "Normal"))
data_use$Type <- colors

jpeg(file=paste(args[3], "/", args[2], ".Z.jpeg", sep=""), width=2000, height=1000)
p <- ggplot(data_use, aes(x=ID, y=Data, color=Type)) + geom_point()
p + theme(axis.text.x=
		element_text(
			size=10,
			vjust=0.5,
			hjust=0.5,
			angle=90
		)
	) +
	labs(x="Location ID", y=args[2]) +
	scale_y_continuous(breaks=c(-5.0, -4.5, -4.0, -3.5, -3.0, 0, 3.0, 3.5, 4.5, 5.0)) +
	geom_rect(mapping=aes(ymin=-3.0, ymax=3.0, xmin=0, xmax=Inf), alpha=0.01, color=NA)
dev.off()