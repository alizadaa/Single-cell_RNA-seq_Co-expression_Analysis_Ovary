library(dplyr)
library(Seurat)
library(patchwork)
set.seed(1234)

# correlations from the combined re-clustered germline clusters

natcomm <- readRDS(file = "/Users/alizad01/Desktop/single_cell/GSE136162_RAW/combined_datasets.rds")
Cluster.10 <- subset(x = natcomm, idents = c('10'))
cl10counts <- Cluster.10[["RNA"]]@counts
matrix_mod<-as.matrix(cl10counts)
list = c("AGO3","aub","vas","qin","rhi","del","cuff","Nxf3","tej","CG13741","CG12721","spn-E","armi","arx","BoYb","emb","CG10880","fs(1)Yb","Gasz","Hen1","krimp","mael","mino","Nbr","nxf2","CG5694","Panx","papi","piwi","shu","SoYb","squ","tapas","thoc5","thoc7","Trf2","tud","vls","Su(var)2-10","Hel25E","egg","tsu","CG14438","csul")

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "germline_cor_GSE136162_combined.txt", quote =FALSE, sep="\t",col.names=list)

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)$p.value})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "germline_cor_GSE136162_combined_pvals.txt", quote =FALSE, sep="\t",col.names=list)


plosbio <- readRDS(file = "/single_cell/GSE146040_RAW/sep_GSE146040.rds")
Cluster.3 <- subset(x = plosbio, idents = c('3'))
cl3counts <- Cluster.3[["RNA"]]@counts
matrix_mod<-as.matrix(cl3counts)

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "germline_cor_GSE146040.txt", quote =FALSE, sep="\t",col.names=list)

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)$p.value})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "germline_cor_GSE146040_pvals.txt", quote =FALSE, sep="\t",col.names=list)

larva <- readRDS(file = "/single_cell/GSE131971_RAW/combined.rds")
Cluster.0 <- subset(x = natcomm, idents = c('0'))
cl0counts <- Cluster.0[["RNA"]]@counts
matrix_mod<-as.matrix(cl0counts)
list = c("AGO3","aub","vas","qin","rhi","del","cuff","Nxf3","tej","CG13741","CG12721","spn-E","armi","arx","BoYb","emb","CG10880","fs(1)Yb","Gasz","Hen1","krimp","mael","mino","Nbr","nxf2","CG5694","Panx","papi","piwi","shu","SoYb","Nos","squ","tapas","thoc5","thoc7","Trf2","tud","vls","Su(var)2-10","Hel25E","egg","tsu","CG14438","csul")

c=lapply(list, function(x) {
  gene <- as.numeric(matrix_mod[x,])
  correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/mnt/scratcha/ghlab/alizada/single_cell/all_correlations/germline_cor_GSE131971_combined.txt", quote =FALSE, sep="\t",col.names=list)


c=lapply(list, function(x) {
  gene <- as.numeric(matrix_mod[x,])
  correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)$p.value})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/mnt/scratcha/ghlab/alizada/single_cell/all_correlations/germline_cor_GSE131971_combined_pvals.txt", quote =FALSE, sep="\t",col.names=list)

