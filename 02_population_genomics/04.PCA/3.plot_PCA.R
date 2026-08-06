# 加载必要的包
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggforce) # 用于绘制圈

# 读取数据
eigenval <- scan("pca.eigenval")  # 读取特征值文件
eigenvec <- read.table("pca.eigenvec", header = T)  # 读取特征向量文件
sample_info <- read.table("pigeon119_filtered_biAllele_missing1_chrAuto_maf0.01.recode.prune0.1.group_3", 
                          header = FALSE, col.names = c("ID", "Breed","Class"))  # 读取品种信息

# 计算每个主成分的解释百分比
variance_explained <- eigenval / sum(eigenval) * 100

# 准备绘图数据
pca_data <- data.frame(
  Sample = eigenvec[, 2],  # 假设第二列是样本ID
  PC1 = eigenvec[, 3],    # 第三列通常是PC1
  PC2 = eigenvec[, 4],    # 第四列通常是PC2
  Breed = sample_info$Breed[match(eigenvec[, 2], sample_info$ID)],
  Class = factor(sample_info$Class[match(eigenvec[, 2], sample_info$ID)])
)


# 创建25种颜色的调色板
# 方法1: 使用RColorBrewer的Set3调色板(12色) + Set1调色板(9色) + Dark2调色板(8色)
#colors_25 <- c(brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(8, "Dark2")[1:4])
#方法3: 手动指定25种颜色
colors_25 <- c("#BCBD22", "#9370DB", "#8C6D31", "#D62728", "#CC0000",
               "#C71585", "#D2B48C", "#8A2BE2", "#FF6347", "#9467BD",
               "#800080", "#8DC3DB", "#98DF8A", "#5C3317", "#C5B0D5",
               "#FF1493", "#2959A0", "#FF4500", "#DBDB8D", "#9E0505",
               "#4B7CB8", "#637939", "#2CA02C", "#BC8F8F", "#1C2B53")

# 为分类创建形状映射（1=圆形, 2=三角形, 3=正方形）
shape_mapping <- c("1" = 16, "2" = 17, "3" = 15)

# 绘制PCA图  
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Breed, shape = Class)) +
  geom_point(size = 7, alpha = 0.6) +
  labs(
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
    title = "PCA Plot Colored by Breed",
    color = "Breed",
    shape = "Class"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = "black", size = 1.2),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8)  # 缩小图例文字以适应更多品种
  ) +
  scale_color_manual(values = colors_25) +  # 使用自定义的25种颜色
  guides(color = guide_legend(ncol = 2))  # 将图例分为两列显示

# 显示图形
print(p)

# 保存图形
ggsave("PCA_plot_by_breed-20250707.pdf", width = 12, height = 8)
