library(ggplot2)

reads <- read.csv('reads_length.txt', sep=' ', header=FALSE)

ggplot(reads, aes(x=reads$V2, y=reads$V1)) + 
	geom_bar(stat='identity') + 
	xlab('read length') + 
	ylab('counts') + 
	ggtitle('Read Length Distribution')
