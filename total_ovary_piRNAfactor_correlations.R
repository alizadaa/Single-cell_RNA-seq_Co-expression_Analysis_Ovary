library(dplyr)
library(Seurat)
library(patchwork)
set.seed(1234)

pbmc.data1 <- Read10X(data.dir = "/single_cell/GSE136162_RAW/dataset1/")
pbmc.data2 <- Read10X(data.dir = "/single_cell/GSE136162_RAW/dataset2/")
pbmc.data3 <- Read10X(data.dir = "/single_cell/GSE136162_RAW/dataset3/")

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

merged <- cbind.fill(pbmc.data1, pbmc.data2, pbmc.data3)
matrix_mod<-as.matrix(merged)

list = c("AGO3","aub","vas","qin","rhi","del","cuff","Nxf3","tej","CG13741","CG12721","spn-E","armi","arx","BoYb","emb","CG10880","fs(1)Yb","Gasz","Hen1","krimp","mael","mino","Nbr","nxf2","CG5694","Panx","papi","piwi","shu","SoYb","squ","tapas","thoc5","thoc7","Trf2","tud","vls","Su(var)2-10","Hel25E","egg","tsu","CG14438","csul")

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/Adult_ovary1_GSE136162_RAW_piRNAfactors_correlation_values.txt", quote =FALSE, sep="\t",col.names=list)

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)$p.value})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/Adult_ovary1_GSE136162_cor_pvalues.txt", quote =FALSE, sep="\t",col.names=list)

pbmc.data4 <- Read10X(data.dir = "/single_cell/GSE146040/")
matrix_mod<-as.matrix(pbmc.data4 )

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/Adult_ovary2_GSE146040_piRNAfactors_correlation_values.txt", quote =FALSE, sep="\t",col.names=list)

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)$p.value})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/Adult_ovary2_GSE146040_cor_pvalues.txt", quote =FALSE, sep="\t",col.names=list)

pbmc.data5 <- Read10X(data.dir = "/single_cell/GSE131971_RAW/dataset1/")
pbmc.data6 <- Read10X(data.dir = "/single_cell/GSE131971_RAW/dataset2/")

merged <- cbind.fill(pbmc.data5, pbmc.data6)
matrix_mod<-as.matrix(merged)

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/Larva_ovary_GSE131971_piRNAfactors_correlation_values.txt", quote =FALSE, sep="\t",col.names=list)

c=lapply(list, function(x) {
gene <- as.numeric(matrix_mod[x,])
correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)$p.value})
})
b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/Larva_ovary_GSE131971_cor_pvals.txt", quote =FALSE, sep="\t",col.names=list)

