setwd('~/Desktop/514/')
library(qqman)

# Load the rall.txt file into RStudio
SCZ <- read.table('~/Desktop/rall.txt', header = TRUE, sep = "\t")

#Display the structure of SCZ
str(SCZ)

#To see how all chromosomes in SCZ are currently labeled
levels(SCZ$hg19chrc)

#Remove “chr” from chromosome names
SCZ$hg19chrc <- gsub("chr","",SCZ$hg19chrc)

head(SCZ)
tail(SCZ)

#Change chromosome X to chromosome 23
SCZ$hg19chrc <- gsub("X","23",SCZ$hg19chrc)

tail(SCZ)

#Check to see if the chromosome column is set as an integer
is.integer(SCZ$hg19chrc)
#Set chromosome as Integer
SCZ$hg19chrc <- as.integer(SCZ$hg19chrc)
is.integer(SCZ$hg19chrc)

str(SCZ)

#change the first [1] column’s name to “CHR”
colnames(SCZ)[1] <- "CHR"
#change the fifth [5] column’s name to “BP”
colnames(SCZ)[5] <- "BP"
#change the ninth [9] column’s name to “P”
colnames(SCZ)[9] <- "P"

str(SCZ)

#Calculate the number of the significance level
0.05/9444230

manhattan(SCZ, chrlabs = c(1:22, "X"), suggestiveline = FALSE,
          genomewideline = -log10(5.3e-09), main = "PGC SCZ GWAS")

manhattan(subset(SCZ, CHR == 6), suggestiveline = FALSE,
          genomewideline = -log10(5.3e-09), main = "PGC SCZ GWAS chr6")
