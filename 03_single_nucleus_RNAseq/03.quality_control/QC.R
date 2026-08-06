library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(harmony)

setwd("/project-whj/sis/2.tumbler/10.SC/3.Seurat+scrublet/2.qc_2")

# 设置根目录路径
root_path <- '/project-whj/sis/2.tumbler/10.SC/1.expression_matrix'

# 列出根目录下的所有子目录
subdirectories <- list.dirs(path = root_path, recursive = FALSE)

# 使用正则表达式筛选出包含的目录
selected_subdirectories <- subdirectories[grepl(".*WMS.*", subdirectories)]

# 构造完整路径
file_names <- "filtered_cell_gene_matrix"
file_paths <- file.path(selected_subdirectories, "outs", file_names)

#循环读取数据，并生成seurat对象
sceList = lapply(file_paths,function(folder){
  CreateSeuratObject(counts = Read10X(folder),
                     project = folder)
})

#样品信息
cell.ids = basename(selected_subdirectories)

#将所有样品合并成一个seurat对象
pigeon_brain <- merge(sceList[[1]],
                      y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]],sceList[[10]],sceList[[11]],sceList[[12]],sceList[[13]],sceList[[14]],sceList[[15]],sceList[[16]],sceList[[17]],sceList[[18]],sceList[[19]],sceList[[20]],sceList[[21]],sceList[[22]],sceList[[23]],sceList[[24]]),
                      add.cell.ids = cell.ids,
                      project = "pigeon_brain")

table(pigeon_brain@meta.data$orig.ident)
rownames(pigeon_brain)
colnames(pigeon_brain)

#修改orig.ident
meta.data <- pigeon_brain@meta.data
meta.data$orig.ident <- sapply(X = strsplit(meta.data$orig.ident, split = "/"), FUN = "[", 7)
pigeon_brain <- AddMetaData(pigeon_brain, meta.data)


#####增加个体信息和样品种类信息
meta.data <- pigeon_brain@meta.data
meta.data$individual.id = sapply(X = strsplit(meta.data$orig.ident, split = "-"), FUN = "[", 1)
meta.data$brain.region = sapply(X = strsplit(meta.data$orig.ident, split = "-"), FUN = "[", 2)
pigeon_brain <- AddMetaData(pigeon_brain, meta.data)

#数据的基本情况--小提琴图
genes_pigeon <- as.data.frame(rownames(pigeon_brain))
mt<-c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","CYTB","ND6")
found_genes <- mt %in% genes_pigeon$`rownames(pigeon_brain)`
print(found_genes)
pigeon_brain[["percent.mt"]] <- PercentageFeatureSet(pigeon_brain, features=mt) #提取有关线粒体的基因
Idents(pigeon_brain) <- "orig.ident"
ribo<-c("RPL22","RPL11","RPS6KA1","RPS8","RPL5","RPS27","RPS6KC1","RPS7","RPS27A","RPL31","RPL37A","RPL32","RPL15","RPSA","RPL14")
found_genes <- ribo %in% genes_pigeon$`rownames(pigeon_brain)`
print(found_genes)
pigeon_brain[["percent.ribo"]] <- PercentageFeatureSet(pigeon_brain, features=ribo) #提取有关核糖体的基因

# 基因数量随细胞数量的变化关系
plot1 <- FeatureScatter(pigeon_brain, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.4)
plot2 <- FeatureScatter(pigeon_brain, "nCount_RNA", "percent.mt", group.by = "orig.ident", pt.size = 0.4)
plot3 <- FeatureScatter(pigeon_brain, "nCount_RNA", "percent.ribo", group.by = "orig.ident", pt.size = 0.4)
pdf(file = "1.VlnPlot_rawdata_nFeatur_nCount_percent.mt.ribo.pdf",width=10,height=16)
print(VlnPlot(object = pigeon_brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 1,pt.size = 0.4))
print(VlnPlot(object = pigeon_brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 1,pt.size = 0))
grid.arrange(plot1, plot2, plot3, nrow = 3)
dev.off()

#数据的基本情况--直方图
cell_data = as.data.frame(pigeon_brain@meta.data)
pdf(file ="2.Histogram_rawdata_nFeatur_nCount.pdf",width=10,height=8)
ggplot(data=cell_data,aes(x=nFeature_RNA)) +labs(x="nFeature_RNA",y="Counts")+
  geom_histogram(fill='#6495ED',alpha=0.5,bins = 100)+
  scale_x_continuous(breaks = seq(from = 0, to = 9000, by = 500))
ggplot(data=cell_data,aes(x=nCount_RNA)) +labs(x="nCount_RNA",y="Counts")+
  geom_histogram(fill='#FFA500',alpha=0.5,bins = 100)+coord_cartesian(xlim = c(0,90000))
ggplot(data=cell_data,aes(x=percent.mt)) +labs(x="percent.mt",y="Counts")+
  geom_histogram(fill='#FFA500',alpha=0.5,bins = 30)
ggplot(data=cell_data,aes(x=percent.ribo)) +labs(x="percent.ribo",y="Counts")+
  geom_histogram(fill='#FFA500',alpha=0.5,bins = 30)
dev.off()

#####-------------------QC
#####-------------------质控2：去除双胞
sample_ids <- unique(pigeon_brain@meta.data$orig.ident)
scrublet <- list()

# 读取 scrublet 结果
for (sample_id in sample_ids) {
  file_path <- paste0("/project-whj/sis/2.tumbler/10.SC/3.Seurat+scrublet/1.scrublet/", sample_id, "-doublet.txt")
  scrublet_result <- read.table(file_path, header = T, sep = ",", stringsAsFactors = FALSE)
  scrublet_result$barcode <- paste0(sample_id, "_", scrublet_result$barcode)
  scrublet[[sample_id]] <- scrublet_result
}

# 初始化 Seurat 对象的 scrublet 列为 "Singlet"
pigeon_brain@meta.data$scrublet <- "Singlet"

# 更新 Seurat 对象的 scrublet 列
for (sample_id in sample_ids) {
  scrublet_result <- scrublet[[sample_id]]
  sample_cells <- rownames(pigeon_brain@meta.data)[pigeon_brain@meta.data$orig.ident == sample_id]
  scrublet_barcodes <- scrublet_result$barcode
  scrublet_calls <- ifelse(scrublet_result$predicted_doublet == "True", "Doublet", "Singlet")
  pigeon_brain@meta.data[sample_cells, "scrublet"] <- scrublet_calls[match(sample_cells, scrublet_barcodes)]
  # 检查是否有未更新的细胞
  if (any(is.na(pigeon_brain@meta.data$scrublet))) {
    cat("存在未更新的细胞，检查匹配结果\n")
  }
}

# 筛选出单细胞（Singlet）
singlet_cells <- rownames(pigeon_brain@meta.data)[pigeon_brain@meta.data$scrublet == "Singlet"]

# 创建一个新的 Seurat 对象，只包含单细胞
pigeon_qc <- subset(pigeon_brain, cells = singlet_cells)

# 统计过滤信息
filter_info <- as.data.frame(table(pigeon_brain$orig.ident))
colnames(filter_info) <- c("sample", "raw")
filter_single <- as.data.frame(table(pigeon_qc$orig.ident))
colnames(filter_single) <- c("sample", "filter_single")
filter_info$filter_single <- filter_single$filter_single[match(filter_info$sample, filter_single$sample)]

# 保存处理后的Seurat对象
saveRDS(pigeon_qc, file = "./pigeon_qc.rds")


pigeon_qc <- readRDS("../2.qc/pigeon_qc.rds")
#####-------------------质控2：去除低质量细胞和潜在双胞
filter_info <- as.data.frame(table(pigeon_qc$orig.ident))
colnames(filter_info) <- c("sample", "filter_single")

cell_data = as.data.frame(pigeon_qc@meta.data)
pigeon_qc <- subset(x = pigeon_qc, subset = percent.mt <= 5)
filter_mito <- as.data.frame(table(pigeon_qc$orig.ident))
colnames(filter_mito) <- c("sample", "filter_mito")
filter_info$filter_mito <- filter_mito$filter_mito[match(filter_info$sample, filter_mito$sample)]

median_nFeature_DN <- cell_data %>% filter(brain.region == "DN") %>%
  summarise(median = median(nFeature_RNA)) %>% pull(median)
std_deviation_nFeature_DN <- cell_data %>% filter(brain.region == "DN") %>%
  summarise(sd = sd(nFeature_RNA)) %>% pull(sd)
max_nFeature_RNA_DN = median_nFeature_DN + 3*std_deviation_nFeature_DN

median_nFeature_XN <- cell_data %>% filter(brain.region == "XN") %>%
  summarise(median = median(nFeature_RNA)) %>% pull(median)
std_deviation_nFeature_XN <- cell_data %>% filter(brain.region == "XN") %>%
  summarise(sd = sd(nFeature_RNA)) %>% pull(sd)
max_nFeature_RNA_XN = median_nFeature_XN + 3*std_deviation_nFeature_XN

median_nCount_DN <- cell_data %>% filter(brain.region == "DN") %>%
  summarise(median = median(nCount_RNA)) %>% pull(median)
std_deviation_nCount_DN <- cell_data %>% filter(brain.region == "DN") %>%
  summarise(sd = sd(nCount_RNA)) %>% pull(sd)
max_nCount_RNA_DN = median_nCount_DN + 3*std_deviation_nCount_DN

median_nCount_XN <- cell_data %>% filter(brain.region == "XN") %>%
  summarise(median = median(nCount_RNA)) %>% pull(median)
std_deviation_nCount_XN <- cell_data %>% filter(brain.region == "XN") %>%
  summarise(sd = sd(nCount_RNA)) %>% pull(sd)
max_nCount_RNA_XN = median_nCount_XN + 3*std_deviation_nCount_XN


# 过滤条件
filter_brain_region <- function(seurat_object, region_name, nFeature_min, nFeature_max, nCount_min, nCount_max) {
  subset(seurat_object, subset = brain.region == region_name &
           nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max &
           nCount_RNA > nCount_min & nCount_RNA < nCount_max)
}

# 对不同的 brain.region 应用过滤
pigeon_brain_DN <- filter_brain_region(pigeon_qc, "DN", 200, max_nFeature_RNA_DN, 1000, max_nCount_RNA_DN)
pigeon_brain_XN <- filter_brain_region(pigeon_qc, "XN", 200, max_nFeature_RNA_XN, 1000, max_nCount_RNA_XN)

# 合并过滤后的结果
filtered_pigeon_brain <- merge(
  x = pigeon_brain_DN,
  y = pigeon_brain_XN,
  add.cell.ids = c("DN", "XN")
)
filter_qc <- as.data.frame(table(filtered_pigeon_brain$orig.ident))
colnames(filter_qc) <- c("sample", "filter_qc")
filter_info$filter_qc <- filter_qc$filter_qc[match(filter_info$sample, filter_qc$sample)]

#质控后再次进行作图
plot1 <- VlnPlot(filtered_pigeon_brain, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
plot2 <- VlnPlot(filtered_pigeon_brain, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
plot3 <- VlnPlot(filtered_pigeon_brain, features = "percent.mt", pt.size = 0.1) + NoLegend()
plot1.1 <- VlnPlot(filtered_pigeon_brain, features = "nFeature_RNA", pt.size = 0) + NoLegend()
plot2.1 <- VlnPlot(filtered_pigeon_brain, features = "nCount_RNA", pt.size = 0) + NoLegend()
plot3.1 <- VlnPlot(filtered_pigeon_brain, features = "percent.mt", pt.size = 0) + NoLegend()
plot1.2 <- VlnPlot(pigeon_brain_DN, features = "nFeature_RNA", pt.size = 0) + NoLegend()
plot2.2 <- VlnPlot(pigeon_brain_XN, features = "nFeature_RNA", pt.size = 0) + NoLegend()
plot1.3 <- VlnPlot(pigeon_brain_DN, features = "nCount_RNA", pt.size = 0) + NoLegend()
plot2.3 <- VlnPlot(pigeon_brain_XN, features = "nCount_RNA", pt.size = 0) + NoLegend()
plot1.4 <- VlnPlot(pigeon_brain_DN, features = "percent.mt", pt.size = 0) + NoLegend()
plot2.4 <- VlnPlot(pigeon_brain_XN, features = "percent.mt", pt.size = 0) + NoLegend()
pdf(file = "3.violin_plots_filtered.pdf", width = 10, height = 16)
grid.arrange(plot1, plot2, plot3, nrow = 3)
grid.arrange(plot1.1, plot2.1, plot3.1, nrow = 3)
grid.arrange(plot1.2, plot1.3, plot1.4, nrow = 3)
grid.arrange(plot2.2, plot2.3, plot2.4, nrow = 3)
dev.off()

write.table(filter_info,file = 'filter_info.txt',sep = '\t', quote = F)
saveRDS(filtered_pigeon_brain,file = "filtered_pigeon_brain.rds")
saveRDS(pigeon_brain_DN,file = "filtered_pigeon_brain_DN.rds")
saveRDS(pigeon_brain_HM,file = "filtered_pigeon_brain_HM.rds")
