library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(harmony)

setwd("/project-whj/sis/2.tumbler/10.SC/3.Seurat+scrublet/3.seurat_DN/fi")

#####-------------------标准化-归一化-高变基因（耗内存）
scRNA<-readRDS("/project-whj/sis/2.tumbler/10.SC/3.Seurat+scrublet/2.qc_2/filtered_pigeon_brain_DN.rds")
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000, verbose=FALSE)
##高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000, verbose=FALSE)
P1 <- VariableFeaturePlot(scRNA)
P1 <- P1 +   theme( text = element_text(size = 12))
top10_genes <- head(VariableFeatures(scRNA), 10)
P2 <- LabelPoints(plot = P1, points = top10_genes, repel = TRUE)
P2 <- P2 +   theme( text = element_text(size = 12))
pdf(file="4.variableFeatures.pdf",width = 12, height = 4)
P1 + P2
dev.off()

##归一化
#对基因在细胞间的表达量进行线性变换，使平均值为0，方差为1，去除高表达基因的影响
scRNA_scaled <- ScaleData(scRNA,features = rownames(scRNA), verbose=FALSE)

#去除线粒体的影响
scRNA_scaled <- ScaleData(scRNA, vars.to.regress = "percent.mt")
dim(scRNA_scaled)
scRNA <- scRNA_scaled


##PCA
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA), verbose=FALSE)
pdf(file="5.PCA.pdf",width = 12, height = 8)
VizDimLoadings(scRNA, dims=1:2, reduction="pca")
DimPlot(scRNA, reduction="pca")
dev.off()
# 碎石图，确定维度 （推荐）
ElbowPlot(scRNA, reduction="pca", ndims=50) 
ggsave("6.确定维度碎石图.png", width=10, height=10, dpi=600)

table(scRNA@meta.data$individual.id)
scRNA@meta.data$group <- NA
# 根据 individual.id 分配组
scRNA@meta.data$group[scRNA@meta.data$individual.id %in% c("WMS420", "WMS421", "WMS422")] <- "homing"
scRNA@meta.data$group[scRNA@meta.data$individual.id %in% c("WMS508", "WMS509", "WMS510")] <- "tumbler"

###   umap_tsne
scRNA <- RunUMAP(scRNA, dims = 1:50,reduction.name = "umap", repel = TRUE)
scRNA <- RunTSNE(scRNA, dims = 1:50, reduction.name = "tsne", repel = TRUE)
pdf(file="7.umap-tsne.pdf", width = 12, height = 12)
DimPlot(scRNA, reduction="umap", group.by = 'orig.ident', label=TRUE)
DimPlot(scRNA, reduction="umap", group.by = 'brain.region',label=TRUE)
DimPlot(scRNA, reduction="umap", group.by = 'group',label=TRUE)
DimPlot(scRNA, reduction="tsne", group.by = 'orig.ident', label=TRUE)
DimPlot(scRNA, reduction="tsne", group.by = 'brain.region',label=TRUE)
DimPlot(scRNA, reduction="tsne", group.by = 'group',label=TRUE)
dev.off()

##去批次效应
scRNA <- scRNA %>% RunHarmony("individual.id", plot_convergence = TRUE)

###   umap_tsne
scRNA <- RunUMAP(scRNA, reduction="harmony", dims = 1:50, )
scRNA <- RunTSNE(scRNA, reduction="harmony", dims = 1:50)
pdf(file="7.umap-tsne-harmony.pdf", width = 12, height = 12)
DimPlot(scRNA, reduction="umap", group.by = 'orig.ident', label=TRUE)
DimPlot(scRNA, reduction="umap", group.by = 'brain.region',label=TRUE)
DimPlot(scRNA, reduction="umap", group.by = 'group',label=TRUE)
DimPlot(scRNA, reduction="tsne", group.by = 'orig.ident', label=TRUE)
DimPlot(scRNA, reduction="tsne", group.by = 'brain.region',label=TRUE)
DimPlot(scRNA, reduction="tsne", group.by = 'group',label=TRUE)
dev.off()


#分群、细胞类型注释
scRNA <- FindNeighbors(scRNA, dims=1:50, verbose=FALSE)
scRNA <- FindClusters(scRNA, resolution=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), verbose=FALSE)
library(clustree)
clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
ggsave("8.cluster图.png", width=16, height=10, dpi=600)

resolutions <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
for (res in resolutions) {
  Idents(scRNA) <- scRNA[[paste0("RNA_snn_res.", res)]]
  pdf(file = paste0("9.UMAP-", res, ".pdf"), width = 12, height = 12)
  print(DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle(paste0("Resolution: ", res)))
  print(DimPlot(scRNA, reduction = "umap", group.by = 'orig.ident', label = TRUE, pt.size = 0.5))
  print(DimPlot(scRNA, reduction = "umap", group.by = 'brain.region', label = TRUE, pt.size = 0.5))
  print(DimPlot(scRNA, reduction = "umap", group.by = 'group', label = TRUE, pt.size = 0.5))
  dev.off()
}

saveRDS(scRNA, file = "scRNA.rds")

scRNA <- readRDS("./scRNA.rds")
# 找到每一个cluster当中的marker，并且只展示阳性的marker。
scRNA <- FindClusters(scRNA, resolution=4, verbose=FALSE)
pigeon.markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, test.use = "wilcox",
                                 min.pct = 0.25, only.pos = TRUE)
pigeon.markers <- pigeon.markers %>% group_by(cluster)
write.table(pigeon.markers,file = 'pigeon.markers.res4.txt',sep = '\t', quote = F)

top5 <- pigeon.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5 <- top5[!duplicated(top5$gene), ]
top5 <- na.omit(top5)
clusters <- unique(pigeon.markers$cluster)

pdf(file = "Cluster_DotPlot_res4.pdf",width=64, height=24)
DotPlot(scRNA, assay='RNA', features = top5$gene) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
saveRDS(scRNA, file = "scRNA_marker_res4.rds")
