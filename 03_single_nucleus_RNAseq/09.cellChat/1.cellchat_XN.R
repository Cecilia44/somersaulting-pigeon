# 或者使用 reticulate 包（如果你在 RStudio 中调用 Python）
library(reticulate)
use_condaenv("sc_RNA", required = TRUE)
use_condaenv("/project-whj/sis/software/anaconda3/envs/sc_RNA", required = TRUE)

library(CellChat)
library(patchwork)
library(Seurat)
library(zellkonverter)

#ad <- import("anndata")
#adata <- ad$read_h5ad("DN_celltype_res4.h5ad")
## 2. 转换为 Seurat 对象
## 注意：AnnData 的 X 矩阵通常对应 Seurat 的 counts 或 data
## 3. 提取矩阵、细胞名和基因名
## 注意：Python 矩阵是 (cell, gene)，R 需要 (gene, cell)，所以要转置 t()
#counts <- t(adata$X) 
#colnames(counts) <- adata$obs_names$to_list()
#rownames(counts) <- adata$var_names$to_list()
#
## 4. 提取元数据 (Metadata)
#metadata <- adata$obs
## 5. 创建并保存 RDS
#scRNA <- CreateSeuratObject(counts = counts, meta.data = metadata)
## 检查一下基因名是否是大写 Symbol
#head(rownames(scRNA))
#
#saveRDS(scRNA, "DN_celltype_res4.rds")

scRNA <- readRDS("../scRNA_marker_res5_remove50_markerAnnotation.rds")

# 1. 确认细胞数量和分群数量
table(Idents(scRNA)) 
# 1. 切换 Ident
Idents(scRNA) <- scRNA$cell.cluster
table(Idents(scRNA))
# 定义需要排除的类群
remove_cells <- c("red blood cells", 
                  "fibroblast", 
                  "endothelial cells", "choroid")

# 执行过滤：只保留不在 remove_cells 列表中的细胞
scRNA_clean <- subset(scRNA, idents = remove_cells, invert = TRUE)

# 验证过滤结果：确认那四个类群已经消失
table(Idents(scRNA_clean))

sc_homing <- subset(scRNA_clean, group == "homing")
sc_tumbler <- subset(scRNA_clean, group == "tumbler")

library(future)
# --- 核心配置 ---
# 解决你之前遇到的 4.18 GiB 限制，设置为 20GB 确保全量数据通过
options(future.globals.maxSize = 80 * 1024^3)
# 针对全量数据，建议 workers 不要超过 4 个，否则物理内存会瞬间爆满
plan("multisession", workers = 10)

# --- 定义计算函数 ---
run_full_cellchat <- function(seurat_obj, group_name) {
  message(paste("开始处理组别:", group_name))
  
  # 1. 创建对象
  cc <- createCellChat(object = seurat_obj, group.by = "cell.cluster")
  cc@DB <- CellChatDB.human 
  
  # 2. 预处理
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  gc() # 强制内存回收
  
  # 3. 推断通讯概率（最耗时的一步，全量数据预计 40-60 分钟）
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  
  # 4. 保存中间结果，防止后续步骤出错
  saveRDS(cc, paste0("cellchat_", group_name, "_full_results.rds"))
  return(cc)
}

# --- 执行计算 ---

# 1. 计算信鸽 (Homing)
cc_homing <- run_full_cellchat(sc_homing, "homing")
rm(sc_homing); gc() # 跑完立刻删除大对象，释放内存给下一组

# 2. 计算翻翻 (Tumbler)
cc_tumbler <- run_full_cellchat(sc_tumbler, "tumbler")
rm(sc_tumbler); gc()
