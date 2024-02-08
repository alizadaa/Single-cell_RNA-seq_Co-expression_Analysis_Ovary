# Single-cell_RNAseq_Correlation_Analysis
Correlations of genes with piRNA pathway genes using single-cell RNA-seq datasets from the Drosophila ovary. Co-expression analysis encompassing total ovary and germline ovarian cluster-specific correlation of transcription factors with the germline piRNA pathway genes

[![DOI](https://zenodo.org/badge/{754750835}.svg)](https://zenodo.org/badge/latestdoi/{754750835})

# Step 1) Total ovary (larva) single-cell RNA-seq data correlations

library(dplyr)
library(Seurat)
library(patchwork)
set.seed(1234)

data5 <- Read10X(data.dir = "/single_cell/GSE131971_RAW/dataset1/")
data6 <- Read10X(data.dir = "/single_cell/GSE131971_RAW/dataset2/")

merged <- cbind.fill(data5, data6)
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

# Step 2) Extract germ cell clusters from the total ovary (larva) single-cell RNA-seq data

library(monocle3)
library(SeuratWrappers)

data5 <- Read10X(data.dir = "/single_cell/GSE131971_RAW/dataset1/")
data5_obj <- CreateSeuratObject(counts = pbmc.data5, project = "data5", min.cells = 5, min.features = 200)

data6 <- Read10X(data.dir = "/single_cell/GSE131971_RAW/dataset2/")
data6_obj <- CreateSeuratObject(counts = pbmc.data6, project = "data6", min.cells = 5, min.features = 200)
merged_obj <- merge(data5_obj, y = c(data6_obj), add.cell.ids = c("data5", "data6"), project = "merged_Data")

VlnPlot(merged_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
merged_obj <- subset(merged_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
merged_obj <- NormalizeData(merged_obj, normalization.method = "LogNormalize", scale.factor = 10000)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(merged_obj), 10)

all.genes <- rownames(merged_obj)
merged_obj <- ScaleData(merged_obj, features = all.genes)
merged_obj <- RunPCA(merged_obj, features = VariableFeatures(object = merged_obj))
print(merged_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(merged_obj, dims = 1:2, reduction = "pca")

merged_obj <- JackStraw(pbmc, num.replicate = 100)
JackStrawPlot(merged_obj, dims = 1:20)
ElbowPlot(merged_obj)

merged_obj <- FindNeighbors(merged_obj, dims = 1:10)
merged_obj <- FindClusters(merged_obj, resolution = 0.75)

merged_obj <- RunUMAP(pbmc, dims = 1:10)
DimPlot(merged_obj, reduction = "umap")

VlnPlot(pbmc, features = c("vas","AGO3","aub"))
FeaturePlot(pbmc, features = c("vas", "AGO3", "aub"))

cl0counts <- Cluster.0[["RNA"]]@counts
saveRDS(pbmc, file = "/single_cell/GSE136162_RAW/dataset2_res075_dim10.rds")

# Find germline Biomarkers
Cluster.0 <- subset(x =merged_obj, idents = c('0'))
cluster0.markers <- FindMarkers(merged_obj, ident.1 = 0, ident.2 = c(0, 3), min.pct = 0.25) 
germ.markers <- FindAllMarkers(merged_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
germ.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Step 3) Re-cluster germ cell cluster

Cluster.0 <- subset(x =merged_obj, idents = c('0'))
Cluster.0 <- RunPCA(Cluster.0], verbose = FALSE)
Cluster.0 <- FindNeighbors(Cluster.0, dims = 1:3)
Cluster.0 <- FindClusters(Cluster.0, resolution = 0.1)
Cluster.0 <- RunUMAP(Cluster.0, dims = 1:3)

new.cluster.ids <- c("undiff_germ1", "diff_germ","undiff_germ2")
names(new.cluster.ids) <- levels(Cluster.0)
Cluster.0 <- RenameIdents(Cluster.0, new.cluster.ids)
DimPlot(Cluster.0, reduction = "umap")

cds <- as.cell_data_set(Cluster.0) 
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 
cds <- learn_graph(cds, use_partition = TRUE)

undiff_germ1 <- subset(x =Cluster.0, idents = c('undiff_germ1'))
cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells( cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = F,cell_size=1.5)

cds <- estimate_size_factors(cds) 
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Cluster.0)
plot_cells( cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = F,cell_size=1.5,genes=c("bam"))
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("bam", "osk", "vas","bgcn")))
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset)

my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("AGO3","aub","vas","qin","rhi","del","cuff","Nxf3","tej","CG13741","piwi","Panx","CG12721","krimp","zuc")))
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, nrow=3,ncol=5,vertical_jitter=T)

saveRDS(Cluster.0,file = "/single_cell/GSE131971_RAW/germline_cluster_LL3.rds"

# Step 4) Correlations using the combined and re-clustered germline cluster

correlations from the combined re-clustered germline clusters

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
write.table(t, file = "/single_cell/GSE131971_RAW/germline_cor_GSE131971_combined.txt", quote =FALSE, sep="\t",col.names=list)


c=lapply(list, function(x) {
  gene <- as.numeric(matrix_mod[x,])
  correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)$p.value})
})

b = c %>% bind_rows()

t=t(b)
write.table(t, file = "/single_cell/GSE131971_RAW/germline_cor_GSE131971_combined_pvals.txt", quote =FALSE, sep="\t",col.names=list)


