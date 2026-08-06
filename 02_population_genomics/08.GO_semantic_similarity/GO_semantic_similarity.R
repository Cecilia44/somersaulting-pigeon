# -----------------------------------------------------------------------
# 功能：5 组基因集全 GO 集合语义相似度热图 (不区分类别)
# -----------------------------------------------------------------------
library(GOSemSim)
library(org.Gg.eg.db)
library(dplyr)
library(pheatmap)

# 1. 初始化语义数据库 (针对家鸡 BP 过程)
message(">> 正在初始化 GO 语义环境...")
d <- godata('org.Gg.eg.db', ont="BP", computeIC=TRUE, annoDb='org.Gg.eg.db')

# 2. 定义文件列表并提取显著 GO ID
# 请确保文件名与您的本地文件一致
file_list <- list(
  "Common" = "gprofiler_common_11_GOBP.csv",
  "Tumbler_Specific" = "gprofiler_xinjiang_specific140_GOBP.csv",
  "Parlor_Specific" = "gprofiler_parlorRoller_specific74_GOBP.csv",
  "Tumbler_Total" = "gprofiler_Xinjiang_total151_GOBP.csv",
  "Parlor_Total" = "gprofiler_ParlorRoller_total85_GOBP.csv"
)

# 提取各组所有显著 (P.adj < 0.05) 的 GO ID
go_sets <- list()
for(name in names(file_list)) {
  if(file.exists(file_list[[name]])) {
    ids <- read.csv(file_list[[name]]) %>%
      filter(adjusted_p_value < 0.05) %>%
      pull(term_id)
    
    if(length(ids) > 0) {
      go_sets[[name]] <- unique(ids)
      message(paste(">>", name, "提取到", length(ids), "个显著 GO 条目"))
    } else {
      message(paste(">> 警告:", name, "没有显著富集条目"))
    }
  }
}

# 3. 构建 5x5 相似度矩阵
set_names <- names(go_sets)
n <- length(set_names)
sim_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(set_names, set_names))

message(">> 正在计算组间全局语义相似度 (Wang's Method)...")
for(i in 1:n) {
  for(j in i:n) {
    # 计算集合 i 和 集合 j 之间的相似度
    # measure="Wang" 基于图结构，combine="BMA" 计算集合间平均最佳匹配得分
    res <- mgoSim(go_sets[[set_names[i]]], go_sets[[set_names[j]]], 
                  semData = d, measure = "Wang", combine = "BMA")
    
    sim_mat[i,j] <- res
    sim_mat[j,i] <- res
  }
}

# 4. 绘制热图
# ---------------------------------------------------------
pdf("Heatmap_Global_GO_Similarity.pdf", width = 8, height = 7)
pheatmap(sim_mat, 
         main = "Global Functional Similarity between Gene Sets",
         color = colorRampPalette(c("white", "#F39B7FFF", "#DC0000FF"))(100),
         display_numbers = TRUE,     # 在方格内显示相似度数值
         number_format = "%.3f",     # 保留三位小数
         fontsize_number = 10,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         border_color = "white",
         angle_col = 45)             # X 轴标签倾斜 45 度以便阅读
dev.off()

message(">> 分析完成！热图已保存为 Heatmap_Global_GO_Similarity.pdf")
