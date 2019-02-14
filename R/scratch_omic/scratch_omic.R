# install.packages("rvest")
library("rvest")
library("stringr")

url = "https://omictools.com/single-cell-rna-seq-category"
web <- read_html(url, silent=TRUE)

h2 <- web %>% html_nodes("h2") %>% html_text()
h3 <- web %>% html_nodes("h3") %>% html_text()

h2 = trimws(h2)
h3 = trimws(h3)

write.table(h2, "h2.txt", row.names=F, col.names=F, quote=F)
write.table(h3, "h3.txt", row.names=F, col.names=F, quote=F)