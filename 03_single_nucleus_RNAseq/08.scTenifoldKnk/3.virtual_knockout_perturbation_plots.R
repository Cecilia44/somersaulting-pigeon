suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(tibble)
})

# =========================
# 0. 用户参数
# =========================
golgi_file <- "GLRX2_virtualKO_homing_parallel/cerebellum_homing__golgi__diffRegulation.csv"
purkinje_file <- "GLRX2_virtualKO_homing_parallel/cerebellum_homing__purkinje_cells__diffRegulation.csv"

outdir <- "vk_perturbation_figures_positiveZ"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

top_n <- 20
padj_cutoff <- 0.05
use_padj_filter <- TRUE      # 若表中有 padj/FDR 等列，则启用
require_positive_Z <- TRUE   # 只保留 Z > 0
min_distance <- 0            # 可改成更严格阈值，例如 1e-4

# =========================
# 1. 辅助函数
# =========================
find_existing_file <- function(path) {
  candidates <- c(
    path,
    file.path(".", path),
    file.path("data", path),
    file.path("input", path),
    file.path("results", path)
  )
  candidates <- unique(candidates[candidates != ""])
  idx <- file.exists(candidates)
  if (any(idx)) return(normalizePath(candidates[idx][1], winslash = "/"))
  
  # 递归搜索
  base <- basename(path)
  hits <- list.files(".", pattern = paste0("^", gsub("\\.", "\\\\.", base), "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) > 0) return(normalizePath(hits[1], winslash = "/"))
  
  stop("找不到文件: ", path)
}

standardize_columns <- function(df) {
  names(df) <- trimws(names(df))
  low <- tolower(names(df))
  
  gene_col <- names(df)[match(TRUE, low %in% c(
    "gene", "genes", "symbol", "gene_name", "genename", "feature"
  ))]
  
  z_col <- names(df)[match(TRUE, low %in% c(
    "z", "zscore", "z_score", "z.score"
  ))]
  
  dist_col <- names(df)[match(TRUE, low %in% c(
    "distance", "dist", "perturbation", "delta", "score"
  ))]
  
  padj_col <- names(df)[match(TRUE, low %in% c(
    "padj", "p_adj", "p.adj", "fdr", "adj_p_val", "adj_pvalue",
    "adjusted_p_value", "p_val_adj", "qvalue", "q_value"
  ))]
  
  p_col <- names(df)[match(TRUE, low %in% c(
    "pvalue", "p_value", "p.value", "pval", "p_val"
  ))]
  
  sig_col <- names(df)[match(TRUE, low %in% c(
    "significant", "sig", "is_significant"
  ))]
  
  if (is.na(gene_col)) stop("未找到 gene 列，请检查列名。")
  if (is.na(z_col)) stop("未找到 Z 列，请检查列名。")
  if (is.na(dist_col)) stop("未找到 distance 列，请检查列名。")
  
  out <- df %>%
    dplyr::rename(
      gene = !!gene_col,
      Z = !!z_col,
      distance = !!dist_col
    )
  
  if (!is.na(padj_col)) out <- out %>% dplyr::rename(padj = !!padj_col)
  if (!is.na(p_col))    out <- out %>% dplyr::rename(pvalue = !!p_col)
  if (!is.na(sig_col))  out <- out %>% dplyr::rename(significant = !!sig_col)
  
  out
}

read_diffreg <- function(file, cell_type_name) {
  file2 <- find_existing_file(file)
  message("Reading: ", file2)
  
  df <- data.table::fread(file2, data.table = FALSE)
  df <- standardize_columns(df)
  
  df <- df %>%
    dplyr::mutate(
      gene = as.character(gene),
      gene = trimws(gene),
      gene = toupper(gene),
      Z = suppressWarnings(as.numeric(Z)),
      distance = suppressWarnings(as.numeric(distance)),
      cell_type = cell_type_name
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(Z), !is.na(distance))
  
  if ("padj" %in% names(df)) {
    df <- df %>% dplyr::mutate(padj = suppressWarnings(as.numeric(padj)))
  }
  if ("pvalue" %in% names(df)) {
    df <- df %>% dplyr::mutate(pvalue = suppressWarnings(as.numeric(pvalue)))
  }
  if ("significant" %in% names(df)) {
    if (is.logical(df$significant)) {
      # 保持不变
    } else {
      sig_chr <- tolower(as.character(df$significant))
      df$significant <- sig_chr %in% c("true", "t", "1", "yes", "y")
    }
  }
  
  df %>% dplyr::distinct(gene, .keep_all = TRUE)
}

filter_positive_hits <- function(df,
                                 require_positive_Z = TRUE,
                                 use_padj_filter = TRUE,
                                 padj_cutoff = 0.05,
                                 min_distance = 0) {
  out <- df
  
  if (require_positive_Z) {
    out <- out %>% dplyr::filter(Z > 0)
  }
  
  out <- out %>% dplyr::filter(distance > min_distance)
  
  if (use_padj_filter && "padj" %in% names(out)) {
    out <- out %>% dplyr::filter(!is.na(padj), padj < padj_cutoff)
  } else if (use_padj_filter && "significant" %in% names(out)) {
    out <- out %>% dplyr::filter(significant %in% TRUE)
  }
  
  out %>%
    dplyr::arrange(dplyr::desc(Z), dplyr::desc(distance))
}

get_top_positive <- function(df, n_top = 20) {
  df %>%
    dplyr::arrange(dplyr::desc(Z), dplyr::desc(distance)) %>%
    dplyr::slice_head(n = n_top)
}

save_plot2 <- function(p, filename, width = 7, height = 6) {
  ggplot2::ggsave(
    filename = file.path(outdir, paste0(filename, ".pdf")),
    plot = p, width = width, height = height, units = "in", device = cairo_pdf
  )
  ggplot2::ggsave(
    filename = file.path(outdir, paste0(filename, ".png")),
    plot = p, width = width, height = height, units = "in", dpi = 400
  )
}

# =========================
# 2. 读取数据
# =========================
golgi_raw <- read_diffreg(golgi_file, "Golgi")
purk_raw  <- read_diffreg(purkinje_file, "Purkinje")

golgi <- filter_positive_hits(
  golgi_raw,
  require_positive_Z = require_positive_Z,
  use_padj_filter = use_padj_filter,
  padj_cutoff = padj_cutoff,
  min_distance = min_distance
)

purk <- filter_positive_hits(
  purk_raw,
  require_positive_Z = require_positive_Z,
  use_padj_filter = use_padj_filter,
  padj_cutoff = padj_cutoff,
  min_distance = min_distance
)

all_df <- bind_rows(golgi, purk)

write.csv(golgi, file.path(outdir, "Golgi_positive_hits_filtered.csv"), row.names = FALSE)
write.csv(purk,  file.path(outdir, "Purkinje_positive_hits_filtered.csv"), row.names = FALSE)

# =========================
# 3. Top positive perturbed genes 图
# =========================
plot_top_positive <- function(df, cell_type_name, n_top = 20, point_fill = "#C44E52") {
  top_df <- get_top_positive(df, n_top = n_top) %>%
    dplyr::arrange(Z) %>%
    dplyr::mutate(gene = factor(gene, levels = gene))
  
  p <- ggplot(top_df, aes(x = Z, y = gene)) +
    geom_segment(aes(x = 0, xend = Z, y = gene, yend = gene),
                 linewidth = 0.7, colour = "grey65") +
    geom_point(aes(size = distance),
               shape = 21, fill = point_fill, colour = "black", stroke = 0.25) +
    labs(
      title = paste0(cell_type_name, " top positive perturbed genes"),
      x = "Positive perturbation Z score",
      y = NULL,
      size = "Distance"
    ) +
    scale_size_continuous(range = c(3.5, 10)) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(face = "italic")
    )
  
  list(plot = p, table = top_df)
}

golgi_top <- plot_top_positive(golgi, "Golgi", n_top = top_n, point_fill = "#C44E52")
purk_top  <- plot_top_positive(purk,  "Purkinje", n_top = top_n, point_fill = "#4C72B0")

save_plot2(golgi_top$plot, "Fig_Golgi_top_positive_perturbed_genes", width = 7.2, height = 6.6)
save_plot2(purk_top$plot,  "Fig_Purkinje_top_positive_perturbed_genes", width = 7.2, height = 6.6)

write.csv(golgi_top$table, file.path(outdir, "Golgi_top_positive_perturbed_genes.csv"), row.names = FALSE)
write.csv(purk_top$table,  file.path(outdir, "Purkinje_top_positive_perturbed_genes.csv"), row.names = FALSE)

# =========================
# 4. 手动 curated gene panel
# =========================
panel_genes <- list(
  Mitochondria_OXPHOS = c(
    "COX1", "COX2", "COX3", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "CYTB", "ATP6", "ATP8", "NDUFA1", "NDUFA2", "NDUFA6", "NDUFB1", "NDUFB6",
    "UQCR10", "UQCR11", "COX4I1", "COX5A", "COX6A1", "ATP5F1A", "ATP5F1B"
  ),
  Excitatory_Synapse = c(
    "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRID2", "GRIN1", "GRIN2A", "GRIN2B",
    "DLG2", "DLG4", "NRXN1", "NRXN2", "NRXN3", "NBEA", "CNTNAP2", "CLSTN1", "CLSTN2"
  ),
  Inhibitory_Synapse = c(
    "GABRA1", "GABRA6", "GABRB1", "GABRB2", "GABRG2", "SLC6A1", "GAD1", "GAD2"
  ),
  Signaling_AxonGuidance = c(
    "NTRK2", "ROBO2", "SLIT1", "SLIT2", "SLIT3", "PTPRO", "RPS6KA2", "NRK", "ELMO1"
  ),
  Vesicle_Trafficking_Cytoskeleton = c(
    "EXOC2", "RABGAP1L", "TSNARE1", "KIF5B", "MAP1A", "FIGN", "CLSTN1", "CLSTN2", "NEBL"
  ),
  RNA_PostTranscriptional = c(
    "CNOT2", "RBMS3", "TRA2A", "TENT5A", "UCHL3", "CELF2", "PRMT8"
  ),
  Calcium_Signaling = c(
    "CALB1", "CALM1", "CALM2", "CALM3"
  )
)

panel_df <- tibble::enframe(panel_genes, name = "category", value = "gene") %>%
  tidyr::unnest(gene) %>%
  dplyr::mutate(gene = toupper(gene)) %>%
  dplyr::distinct(category, gene)

# =========================
# 5. Focused functional heatmap
# =========================
panel_hits <- all_df %>%
  dplyr::inner_join(panel_df, by = "gene") %>%
  dplyr::select(cell_type, gene, category, Z, distance, dplyr::everything())

if (nrow(panel_hits) > 0) {
  heat_df <- panel_hits %>%
    dplyr::group_by(category, gene, cell_type) %>%
    dplyr::summarise(Z = max(Z, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(category, gene, cell_type = c("Golgi", "Purkinje"), fill = list(Z = 0)) %>%
    dplyr::mutate(gene_label = paste0(category, " | ", gene))
  
  gene_order <- heat_df %>%
    dplyr::group_by(category, gene, gene_label) %>%
    dplyr::summarise(maxZ = max(Z, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(category, dplyr::desc(maxZ), gene)
  
  heat_df <- heat_df %>%
    dplyr::mutate(
      gene_label = factor(gene_label, levels = rev(unique(gene_order$gene_label))),
      cell_type = factor(cell_type, levels = c("Golgi", "Purkinje")),
      category = factor(category, levels = unique(gene_order$category))
    )
  
  p_heat <- ggplot(heat_df, aes(x = cell_type, y = gene_label, fill = Z)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(
      low = "grey95",
      high = "#B22222",
      name = "Positive Z"
    ) +
    scale_y_discrete(labels = function(x) sub("^.* \\| ", "", x)) +
    labs(
      title = "Focused functional heatmap of positively perturbed genes",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = "grey80"),
      strip.text.y = element_text(face = "bold", angle = 0),
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "italic", size = 9),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  save_plot2(p_heat, "Fig_focused_functional_heatmap_positiveZ", width = 6.2, height = 10.8)
  write.csv(heat_df, file.path(outdir, "Focused_functional_heatmap_positiveZ_table.csv"), row.names = FALSE)
}

# =========================
# 6. Curated gene-panel summary 图
# =========================
panel_summary <- all_df %>%
  dplyr::inner_join(panel_df, by = "gene") %>%
  dplyr::group_by(cell_type, category) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    mean_Z = mean(Z, na.rm = TRUE),
    median_Z = median(Z, na.rm = TRUE),
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(panel_summary) > 0) {
  panel_summary <- panel_summary %>%
    dplyr::mutate(
      category = factor(
        category,
        levels = c(
          "Mitochondria_OXPHOS",
          "Excitatory_Synapse",
          "Inhibitory_Synapse",
          "Signaling_AxonGuidance",
          "Vesicle_Trafficking_Cytoskeleton",
          "Calcium_Signaling",
          "RNA_PostTranscriptional"
        )
      ),
      cell_type = factor(cell_type, levels = c("Golgi", "Purkinje")),
      size_plot = pmin(pmax(n_genes, 1), 10)
    )
  
  p_panel <- ggplot(panel_summary, aes(x = cell_type, y = category)) +
    geom_point(aes(size = size_plot, fill = mean_Z),
               shape = 21, colour = "black", stroke = 0.25) +
    scale_size_continuous(
      range = c(3.5, 12),
      breaks = c(1, 2, 3, 5, 8, 10)
    ) +
    scale_fill_gradient(
      low = "grey95",
      high = "#B22222",
      name = "Mean positive Z"
    ) +
    labs(
      title = "Curated gene-panel summary",
      x = NULL,
      y = NULL,
      size = "Gene count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  save_plot2(p_panel, "Fig_curated_gene_panel_summary_positiveZ", width = 7.0, height = 4.8)
  write.csv(panel_summary, file.path(outdir, "Curated_gene_panel_summary_positiveZ.csv"), row.names = FALSE)
}

# =========================
# 7. 导出 panel 成员详细表
# =========================
panel_members <- all_df %>%
  dplyr::inner_join(panel_df, by = "gene") %>%
  dplyr::arrange(category, cell_type, dplyr::desc(Z), dplyr::desc(distance))

if (nrow(panel_members) > 0) {
  write.csv(panel_members, file.path(outdir, "Curated_panel_member_rankings_positiveZ.csv"), row.names = FALSE)
}

# =========================
# 8. session info 和提示
# =========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

message("Done.")
message("Output directory: ", normalizePath(outdir, winslash = "/"))
