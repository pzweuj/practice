library(ggplot2)
args <- commandArgs(TRUE)


inputfile <- read.csv(args[1], header=TRUE, sep="\t")
data_use <- inputfile[, c(1, 2)]
names(data_use) <- c("ID", "Data")

colors <- ifelse(data_use[2]>=1.2, "Up", ifelse(data_use[2]<=0.8, "Down", "Normal"))
data_use$Type <- colors

jpeg(file=paste(args[3], "/", args[2], ".DQ.jpeg", sep=""), width=2000, height=1000)
p <- ggplot(data_use, aes(x=ID, y=Data, color=Type)) + geom_point()

cols <- c("Up"="blue", "Normal"="green", "Down"="red")
p + theme(axis.text.x=
		element_text(
			size=10,
			vjust=0.5,
			hjust=0.5,
			angle=90
		)
	) +
	labs(x="Location ID", y=args[2]) +
	scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 2.0, 2.6)) +
	geom_rect(mapping=aes(ymin=0.8, ymax=1.2, xmin=0, xmax=Inf), alpha=0.01, color=NA) +
	scale_colour_manual(values=cols)
dev.off()