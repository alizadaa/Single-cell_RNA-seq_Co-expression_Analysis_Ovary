# Single-cell_RNAseq_Correlation_Analysis
Correlations of genes with piRNA pathway genes using single-cell RNA-seq datasets from the Drosophila ovary. Co-expression analysis encompassing total ovary and germline ovarian cluster-specific correlation of transcription factors with the germline piRNA pathway genes

[![DOI](https://zenodo.org/badge/{754750835}.svg)](https://zenodo.org/badge/latestdoi/{754750835})

# 1) Total ovary (larva) single-cell RNA-seq data correlations

library(dplyr)
library(Seurat)
library(patchwork)
set.seed(1234)

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


