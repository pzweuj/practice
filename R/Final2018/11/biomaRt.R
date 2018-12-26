library(biomaRt)

ep <- read.table("~/Desktop/epilepsy.genes.fpkm.txt", header=TRUE)
## ----ensembl1-------------
ensembl=useMart("ensembl")

## ----list datasets-------------
listDatasets(ensembl)

## ----ensembl and choose appropiate dataset-------------
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

##----listFilters shows all filters--------
filter <- listFilters(ensembl)
head(filter)

##----listAttributes shows all attributes--------
attributes <- listAttributes(ensembl)
head(attributes)

##----the getBM function is the main query function in biomaRt----
#----convert Ensembl gene IDs to HUGO Gene Nomenclature Committee gene symbols 
ep_trans <- getBM(attributes=c("ensembl_gene_id",
                               "hgnc_symbol",
                               "chromosome_name",
                               "start_position",
                               "end_position"),
                  filters="ensembl_gene_id",
                  values=ep$gene_id,
                  mart=ensembl)

# merge
colnames(ep)[1] <- "ensembl_gene_id" 
ep_anno <- merge(x=ep, y=ep_trans, by="ensembl_gene_id", all.x=TRUE)

# select
ep_select <- ep_anno[which(ep_anno$chromosome_name == "X" &
                             ep_anno$MTLE1_FPKM > 1 &
                             ep_anno$MTLE2_FPKM > 1 &
                             ep_anno$CTRL3_FPKM > 1), ]

write.table(ep_select, "epilepsy_ChrX_FPKM1.txt",
            quote=FALSE, sep="\t", row.names=FALSE)
