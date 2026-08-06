# 1. 环境准备
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
required_packages <- c("ggplot2", "dplyr", "tidyr", "stringr", "readr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# 2. 读取数据 (确保文件名与你上传的一致)
df_pr <- read_csv("gprofiler_ParlorRoller_total85_GOBP.csv") %>% mutate(Source_Group = "PR_Total")
df_xj <- read_csv("gprofiler_Xinjiang_total151_GOBP.csv") %>% mutate(Source_Group = "XJ_Total")

# 3. 识别共有 ID
shared_ids <- intersect(df_pr$term_id, df_xj$term_id)

# 4. 定义分类逻辑
classify_go_optimized <- function(term_name) {
  term_lower <- str_to_lower(term_name)
  case_when(
    str_detect(term_lower, "neuron|neural|synapse|synaptic|brain|axon|nervous|neurogenesis|neurotransmitter|glutamate|action potential|membrane potential|depolarization|repolarization|behavior|learning|cognition|circadian rhythm|pain|seizure|cerebrum|forebrain|neuroblast|conditioned|synaptic transmission|gliogenesis|schwan|amyloid|neuro|corticospinal|locomotion|motor behavior|rhythmic synaptic") ~ "Neurobiology",
    str_detect(term_lower, "immune|immunity|inflammatory|inflammation|cytokine|interferon|leukocyte|lymphocyte|\\bcomplement\\b|myeloid|neutrophil") ~ "Immunity",
    str_detect(term_lower, "cytoskeleton|actin|microtubule|filament|junction|adhesion|projection organization|assembly|biogenesis") ~ "Cellular Organization",
    str_detect(term_lower, "signaling|signal|pathway|transduction|response to|stimulus|receptor|wnt|chemotaxis|taxis") ~ "Signaling & Response",
    str_detect(term_lower, "development|morphogenesis|differentiation|formation|organogenesis|growth|maturation|regeneration|division|angiogenesis|vasculature|vessel|circulatory|cardiac|heart|muscle|striated|bone|skeletal|epithelial|genitalia|lung|salivary|retina|eye|limb|appendage") ~ "Development & Morphogenesis",
    str_detect(term_lower, "transport|transporter|carrier|channel|\\bion\\b|\\bcation\\b|\\banion\\b|calcium|sodium|potassium|chloride|homeostasis|transmembrane|import|export|secretion|\\bph\\b") ~ "Transport & Homeostasis",
    TRUE ~ "Others"
  )
}

# 5. 指定分类顺序 (按你要求排序)
cat_order <- c("Neurobiology", "Immunity", "Cellular Organization", 
               "Signaling & Response", "Development & Morphogenesis", "Transport & Homeostasis")

# 6. 数据处理
plot_data <- bind_rows(df_pr, df_xj) %>%
  filter(term_id %in% shared_ids) %>%
  mutate(Category = classify_go_optimized(term_name)) %>%
  filter(Category %in% cat_order) %>% 
  mutate(Category = factor(Category, levels = rev(cat_order))) %>% # 因子化保证排序
  mutate(GeneRatio = intersection_size / query_size) %>%
  arrange(Category, negative_log10_of_adjusted_p_value) %>%
  mutate(term_name = factor(term_name, levels = unique(term_name)))

# 7. 计算阴影框背景坐标
bg_rects <- plot_data %>%
  mutate(term_idx = as.numeric(term_name)) %>%
  group_by(Category) %>%
  summarize(
    ymin = min(term_idx) - 0.5,
    ymax = max(term_idx) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(fill_col = rep(c("even", "odd"), length.out = n()))

# ---------------------------------------------------------
# 8. 绘图：策略一双向对比图
# ---------------------------------------------------------
p1 <- ggplot() +
  # 背景阴影
  geom_rect(data = bg_rects, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = fill_col), 
            alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values = c("even" = "grey30", "odd" = "white")) +
  # 引导线
  geom_segment(data = plot_data, aes(x = 0, xend = ifelse(Source_Group == "PR_Total", GeneRatio, -GeneRatio), 
                                     y = term_name, yend = term_name), color = "grey80", linetype = "dotted") +
  # 气泡点
  geom_point(data = plot_data, aes(x = ifelse(Source_Group == "PR_Total", GeneRatio, -GeneRatio), 
                                   y = term_name, size = intersection_size, color = negative_log10_of_adjusted_p_value)) +
  # 中心轴
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  # 侧边分类名称
  geom_text(data = bg_rects, aes(x = -max(plot_data$GeneRatio)*1.25, y = (ymin + ymax)/2, label = Category), 
            angle = 90, size = 3.5, fontface = "bold", color = "grey20") +
  # 配色与坐标轴
  scale_color_gradientn(colors = c("#313695", "#74add1", "#fee090", "#f46d43", "#a50026")) +
  scale_x_continuous(labels = abs, limits = c(-max(plot_data$GeneRatio)*1.4, max(plot_data$GeneRatio)*1.4)) +
  theme_minimal() +
  labs(title = "Mirror Plot: Shared GO Terms (PR vs XJ)",
       subtitle = "Neurobiology Prioritized | Categorized Background Shading",
       x = "Gene Ratio (Left: XJ_Total | Right: PR_Total)", y = NULL, 
       size = "Count", color = "-log10(P-adj)") +
  theme(axis.text.y = element_text(size = 8, color = "black"),
        panel.grid = element_blank())

# 在绘图 theme 中加入这些参数
p1 <- p1 + theme(
    # 增加 Y 轴文字间距
    axis.text.y = element_text(size = 7, lineheight = 0.8),
    # 强化分类标签的视觉效果
    plot.margin = unit(c(1, 1, 1, 3), "cm"), # 增加左侧留白给大类标签
    # 调整图例位置，释放横向空间
    legend.position = "bottom",
    legend.box = "horizontal"
)

# 动态调整高度：根据条目数量自动计算保存高度
n_terms <- length(unique(plot_data$term_name))
dynamic_height <- n_terms * 0.25 + 2 # 每个条目给 0.25 英寸

# ---------------------------------------------------------
# 9. 同时保存为 PNG 和 PDF
# ---------------------------------------------------------
# 保存 PNG (适合快速查看和 PPT)
ggsave("GO_Shared_Final_Categorized.png", plot = p1, width = 12, height = 10, dpi = 300, bg = "white")

# 保存 PDF (矢量图，适合投稿和高清晰排版)
ggsave("GO_Shared_Final_Categorized.pdf", plot = p1, width = 12, height = 10, device = "pdf", bg = "white")

message("图像已成功保存为 PNG 和 PDF 格式。")
