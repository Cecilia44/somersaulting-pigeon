options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(tidyr)
  library(stringr)
})

## =========================================
## 1. 用户参数
## =========================================
## 改成你当前已有输出文件所在目录
main_dir <- "GLRX2_virtualKO_homing_parallel"

## 新的后续分析输出目录
outdir <- file.path(main_dir, "downstream_resumed")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## 通路库物种
pathway_species <- "Homo sapiens"

## 显著性阈值
padj_cutoff <- 0.05

## diffRegulation 中用于排序的列
gene_col <- "gene"
rank_col <- "distance"

## 重点细胞类型
focus_tissue <- "cerebellum_homing"
focus_celltypes <- c("golgi", "purkinje cells")

## 热图每个重点细胞类型保留的重点通路数
top_n_pathways_per_celltype <- 12

## =========================================
## 2. 工具函数
## =========================================
sanitize_filename <- function(x) {
  x <- gsub("[/ ]", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  x
}

flatten_list_columns <- function(df) {
  is_list_col <- vapply(df, is.list, logical(1))
  if (any(is_list_col)) {
    df[is_list_col] <- lapply(df[is_list_col], function(col) {
      vapply(col, function(x) paste(x, collapse = ";"), character(1))
    })
  }
  df
}

safe_fgsea <- function(diff_df, pathways, gene_col = "gene", rank_col = "distance") {
  stopifnot(gene_col %in% colnames(diff_df))
  stopifnot(rank_col %in% colnames(diff_df))

  df <- diff_df %>%
    filter(!is.na(.data[[gene_col]]), !is.na(.data[[rank_col]])) %>%
    distinct(.data[[gene_col]], .keep_all = TRUE)

  stats <- df[[rank_col]]
  names(stats) <- df[[gene_col]]
  stats <- sort(stats, decreasing = TRUE)

  if (length(stats) < 20) return(NULL)

  fgsea(
    pathways = pathways,
    stats = stats,
    minSize = 10,
    maxSize = 500,
    scoreType = "pos"
  ) %>%
    as_tibble() %>%
    arrange(padj, desc(abs(NES)))
}

plot_top_pathways <- function(fgsea_df, title_txt, file, top_n = 15) {
  df <- fgsea_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(NES))) %>%
    slice_head(n = top_n)

  if (nrow(df) == 0) return(NULL)

  p <- ggplot(df, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "NES", title = title_txt) +
    theme_bw(base_size = 12)

  ggsave(file, p, width = 8, height = 6)
}

extract_key_pathways <- function(fgsea_all) {
  key_patterns <- c(
    "OXIDATIVE_PHOSPHORYLATION",
    "REACTIVE_OXYGEN_SPECIES",
    "ROS",
    "OXIDATIVE_STRESS",
    "APOPTOSIS",
    "CELL_DEATH",
    "HYPOXIA",
    "UNFOLDED_PROTEIN_RESPONSE",
    "STRESS",
    "MITOCHON",
    "RESPIRATORY",
    "ATP",
    "ELECTRON_TRANSPORT",
    "SYNAPT",
    "CALCIUM",
    "ION"
  )

  fgsea_all %>%
    filter(grepl(paste(key_patterns, collapse = "|"), pathway, ignore.case = TRUE)) %>%
    arrange(database, padj, desc(abs(NES)))
}

pretty_pathway <- function(x) {
  x %>%
    str_replace_all("^HALLMARK_", "") %>%
    str_replace_all("^GOBP_", "") %>%
    str_replace_all("^REACTOME_", "") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

## =========================================
## 3. 读取已有主结果文件
## =========================================
summary_file <- file.path(main_dir, "GLRX2_virtualKO_homing_celltype_summary_raw.csv")
damage_file <- file.path(main_dir, "GLRX2_virtualKO_homing_damage_ranking.csv")
top_gene_file <- file.path(main_dir, "GLRX2_virtualKO_homing_top20_perturbed_genes_by_celltype.csv")

stopifnot(file.exists(summary_file))
stopifnot(file.exists(damage_file))
stopifnot(file.exists(top_gene_file))

summary_tbl <- read.csv(summary_file, check.names = FALSE)
damage_tbl <- read.csv(damage_file, check.names = FALSE)
top_gene_tbl <- read.csv(top_gene_file, check.names = FALSE)

## =========================================
## 4. 准备通路库
## =========================================
hallmark_df <- msigdbr(species = pathway_species, collection = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

gobp_df <- msigdbr(
  species = pathway_species,
  collection = "C5",
  subcollection = "GO:BP"
)
gobp_list <- split(gobp_df$gene_symbol, gobp_df$gs_name)

reactome_df <- msigdbr(
  species = pathway_species,
  collection = "C2",
  subcollection = "CP:REACTOME"
)
reactome_list <- split(reactome_df$gene_symbol, reactome_df$gs_name)

## =========================================
## 5. 读取已有 diffRegulation 文件并跑 FGSEA
## =========================================
diff_files <- list.files(
  main_dir,
  pattern = "__diffRegulation\\.csv$",
  full.names = TRUE
)

if (length(diff_files) == 0) {
  stop("No per-celltype diffRegulation files found in main_dir.")
}

cat("Found", length(diff_files), "diffRegulation files.\n")

all_fgsea <- list()
failed_files <- character(0)

for (f in diff_files) {
  message("Processing: ", basename(f))

  diff_df <- tryCatch({
    read.csv(f, check.names = FALSE)
  }, error = function(e) {
    message("  Failed to read file: ", e$message)
    failed_files <<- c(failed_files, basename(f))
    return(NULL)
  })

  if (is.null(diff_df)) next

  required_cols <- c("tissue", "cell_type", gene_col, rank_col)
  if (!all(required_cols %in% colnames(diff_df))) {
    message("  Skip: required columns missing in ", basename(f))
    failed_files <- c(failed_files, basename(f))
    next
  }

  tissue_name <- unique(diff_df$tissue)[1]
  ct_name <- unique(diff_df$cell_type)[1]
  base_stub <- paste0(
    sanitize_filename(tissue_name), "__",
    sanitize_filename(ct_name)
  )

  ## ---------- Hallmark ----------
  fg_h <- tryCatch({
    safe_fgsea(diff_df, hallmark_list, gene_col = gene_col, rank_col = rank_col)
  }, error = function(e) {
    message("  Hallmark fgsea failed: ", e$message)
    return(NULL)
  })

  if (!is.null(fg_h)) {
    fg_h <- fg_h %>%
      mutate(
        tissue = tissue_name,
        cell_type = ct_name,
        database = "Hallmark"
      )

    fg_h_export <- flatten_list_columns(fg_h)

    write.csv(
      fg_h_export,
      file.path(outdir, paste0(base_stub, "__Hallmark_fgsea.csv")),
      row.names = FALSE
    )

    plot_top_pathways(
      fg_h,
      paste0("Hallmark | ", tissue_name, " | ", ct_name),
      file.path(outdir, paste0(base_stub, "__Hallmark_top15.pdf")),
      top_n = 15
    )

    all_fgsea[[paste0(base_stub, "__Hallmark")]] <- fg_h_export
  }

  ## ---------- GO BP ----------
  fg_go <- tryCatch({
    safe_fgsea(diff_df, gobp_list, gene_col = gene_col, rank_col = rank_col)
  }, error = function(e) {
    message("  GO BP fgsea failed: ", e$message)
    return(NULL)
  })

  if (!is.null(fg_go)) {
    fg_go <- fg_go %>%
      mutate(
        tissue = tissue_name,
        cell_type = ct_name,
        database = "GO_BP"
      )

    fg_go_export <- flatten_list_columns(fg_go)

    write.csv(
      fg_go_export,
      file.path(outdir, paste0(base_stub, "__GO_BP_fgsea.csv")),
      row.names = FALSE
    )

    plot_top_pathways(
      fg_go,
      paste0("GO BP | ", tissue_name, " | ", ct_name),
      file.path(outdir, paste0(base_stub, "__GO_BP_top15.pdf")),
      top_n = 15
    )

    all_fgsea[[paste0(base_stub, "__GO_BP")]] <- fg_go_export
  }

  ## ---------- Reactome ----------
  fg_re <- tryCatch({
    safe_fgsea(diff_df, reactome_list, gene_col = gene_col, rank_col = rank_col)
  }, error = function(e) {
    message("  Reactome fgsea failed: ", e$message)
    return(NULL)
  })

  if (!is.null(fg_re)) {
    fg_re <- fg_re %>%
      mutate(
        tissue = tissue_name,
        cell_type = ct_name,
        database = "Reactome"
      )

    fg_re_export <- flatten_list_columns(fg_re)

    write.csv(
      fg_re_export,
      file.path(outdir, paste0(base_stub, "__Reactome_fgsea.csv")),
      row.names = FALSE
    )

    plot_top_pathways(
      fg_re,
      paste0("Reactome | ", tissue_name, " | ", ct_name),
      file.path(outdir, paste0(base_stub, "__Reactome_top15.pdf")),
      top_n = 15
    )

    all_fgsea[[paste0(base_stub, "__Reactome")]] <- fg_re_export
  }
}

if (length(all_fgsea) == 0) {
  stop("No FGSEA results generated. Check gene symbols / pathway overlap.")
}

## =========================================
## 6. 汇总 FGSEA 总结果
## =========================================
fgsea_all <- bind_rows(all_fgsea)

write.csv(
  fgsea_all,
  file.path(outdir, "GLRX2_virtualKO_homing_fgsea_all.csv"),
  row.names = FALSE
)

fgsea_sig <- fgsea_all %>%
  filter(!is.na(padj), padj < padj_cutoff) %>%
  arrange(database, padj, desc(abs(NES)))

write.csv(
  fgsea_sig,
  file.path(outdir, "GLRX2_virtualKO_homing_fgsea_significant.csv"),
  row.names = FALSE
)

key_pathway_tbl <- extract_key_pathways(fgsea_all)

write.csv(
  key_pathway_tbl,
  file.path(outdir, "GLRX2_virtualKO_homing_key_pathways_summary.csv"),
  row.names = FALSE
)

top_path_tbl <- fgsea_all %>%
  group_by(tissue, cell_type, database) %>%
  arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top_path_tbl,
  file.path(outdir, "GLRX2_virtualKO_homing_top10_pathways_by_celltype.csv"),
  row.names = FALSE
)

path_count_tbl <- fgsea_all %>%
  mutate(sig = !is.na(padj) & padj < padj_cutoff) %>%
  group_by(tissue, cell_type, database) %>%
  summarise(n_sig_pathway = sum(sig), .groups = "drop")

p_pathcount <- ggplot(
  path_count_tbl,
  aes(x = reorder(paste(tissue, cell_type, database, sep = " | "), n_sig_pathway),
      y = n_sig_pathway)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of significant pathways",
    title = "Pathway impact across cell types after GLRX2 virtual KO"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(outdir, "GLRX2_virtualKO_homing_sig_pathway_counts_barplot.pdf"),
  p_pathcount, width = 11, height = 8
)

## =========================================
## 7. Golgi / Purkinje 专题图
## =========================================

## 7.1 小脑细胞类型损伤图（基于你已有 damage ranking）
damage_plot_df <- damage_tbl %>%
  filter(tissue == focus_tissue) %>%
  mutate(
    highlight = ifelse(cell_type %in% focus_celltypes, "Focus", "Other"),
    label = ifelse(cell_type %in% focus_celltypes, cell_type, "")
  ) %>%
  arrange(damage_score)

p_damage <- ggplot(
  damage_plot_df,
  aes(x = reorder(cell_type, damage_score), y = damage_score, fill = highlight)
) +
  geom_col() +
  coord_flip() +
  geom_text(
    aes(label = label),
    hjust = -0.05,
    size = 4
  ) +
  scale_fill_manual(values = c("Focus" = "#D55E00", "Other" = "grey75")) +
  labs(
    x = NULL,
    y = "Damage score",
    title = "Cell-type vulnerability to GLRX2 virtual knockout in cerebellum"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) +
  expand_limits(y = max(damage_plot_df$damage_score, na.rm = TRUE) * 1.12)

ggsave(
  file.path(outdir, "Figure_damage_ranking_highlight_golgi_purkinje.pdf"),
  p_damage,
  width = 8.5,
  height = 6.5
)

## 7.2 Golgi / Purkinje 重点通路热图
heat_raw <- key_pathway_tbl %>%
  filter(
    tissue == focus_tissue,
    cell_type %in% focus_celltypes
  ) %>%
  group_by(cell_type, pathway) %>%
  summarise(
    NES = NES[which.min(padj)],
    padj = min(padj, na.rm = TRUE),
    database = database[which.min(padj)],
    .groups = "drop"
  ) %>%
  arrange(cell_type, padj, desc(abs(NES)))

top_pathways <- heat_raw %>%
  group_by(cell_type) %>%
  arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = top_n_pathways_per_celltype) %>%
  ungroup() %>%
  pull(pathway) %>%
  unique()

heat_plot_df <- expand.grid(
  cell_type = focus_celltypes,
  pathway = unique(top_pathways),
  stringsAsFactors = FALSE
) %>%
  left_join(
    heat_raw %>% select(cell_type, pathway, NES, padj, database),
    by = c("cell_type", "pathway")
  ) %>%
  mutate(
    pathway_pretty = pretty_pathway(pathway),
    NES = ifelse(is.na(NES), 0, NES),
    cell_type = factor(cell_type, levels = focus_celltypes)
  )

path_order <- heat_plot_df %>%
  group_by(pathway_pretty) %>%
  summarise(max_absNES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_absNES) %>%
  pull(pathway_pretty)

heat_plot_df <- heat_plot_df %>%
  mutate(pathway_pretty = factor(pathway_pretty, levels = path_order))

p_heat <- ggplot(
  heat_plot_df,
  aes(x = cell_type, y = pathway_pretty, fill = NES)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "NES",
    title = "Key enriched pathways in Golgi and Purkinje cells after GLRX2 virtual KO"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9)
  )

ggsave(
  file.path(outdir, "Figure_keypathway_heatmap_golgi_purkinje.pdf"),
  p_heat,
  width = 8.5,
  height = 7.5
)

write.csv(
  heat_plot_df %>% arrange(pathway_pretty, cell_type),
  file.path(outdir, "keypathway_heatmap_data_golgi_purkinje.csv"),
  row.names = FALSE
)

## 7.3 Golgi / Purkinje 各数据库 top pathways
fgsea_focus <- fgsea_sig %>%
  filter(tissue == focus_tissue, cell_type %in% focus_celltypes)

for (ct in focus_celltypes) {
  df_ct <- fgsea_focus %>% filter(cell_type == ct)

  if (nrow(df_ct) == 0) next

  for (db in c("Hallmark", "GO_BP", "Reactome")) {
    df_sub <- df_ct %>%
      filter(database == db) %>%
      arrange(padj, desc(abs(NES)))

    if (nrow(df_sub) == 0) next

    p <- ggplot(
      df_sub %>% slice_head(n = 15),
      aes(x = reorder(pathway, NES), y = NES)
    ) +
      geom_col() +
      coord_flip() +
      labs(
        x = NULL,
        y = "NES",
        title = paste0(ct, " | ", db, " top pathways")
      ) +
      theme_bw(base_size = 12)

    ggsave(
      file.path(outdir, paste0("Figure_", sanitize_filename(ct), "__", db, "_top15.pdf")),
      p,
      width = 8.5,
      height = 6.5
    )
  }
}

## 7.4 Golgi / Purkinje top perturbed genes
top_gene_focus <- top_gene_tbl %>%
  filter(tissue == focus_tissue, cell_type %in% focus_celltypes)

for (ct in focus_celltypes) {
  df_ct <- top_gene_focus %>%
    filter(cell_type == ct) %>%
    arrange(p.adj, desc(distance))

  if (nrow(df_ct) == 0) next

  p_gene <- ggplot(
    df_ct %>% slice_head(n = 20),
    aes(x = reorder(gene, distance), y = distance)
  ) +
    geom_col() +
    coord_flip() +
    labs(
      x = NULL,
      y = "Distance",
      title = paste0("Top perturbed genes in ", ct, " after GLRX2 virtual KO")
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(outdir, paste0("Figure_top_perturbed_genes_", sanitize_filename(ct), ".pdf")),
    p_gene,
    width = 8.5,
    height = 7
  )
}

## =========================================
## 8. 补充导出摘要表
## =========================================
focus_damage_tbl <- damage_tbl %>%
  filter(tissue == focus_tissue, cell_type %in% focus_celltypes) %>%
  arrange(desc(damage_score))

write.csv(
  focus_damage_tbl,
  file.path(outdir, "golgi_purkinje_damage_summary.csv"),
  row.names = FALSE
)

focus_fgsea_summary <- fgsea_focus %>%
  group_by(cell_type, database) %>%
  arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  focus_fgsea_summary,
  file.path(outdir, "golgi_purkinje_top10_pathways_summary.csv"),
  row.names = FALSE
)

writeLines(c(
  paste("main_dir =", main_dir),
  paste("outdir =", outdir),
  paste("pathway_species =", pathway_species),
  paste("padj_cutoff =", padj_cutoff),
  paste("gene_col =", gene_col),
  paste("rank_col =", rank_col),
  paste("focus_tissue =", focus_tissue),
  paste("focus_celltypes =", paste(focus_celltypes, collapse = ", ")),
  paste("n_diff_files =", length(diff_files)),
  paste("n_failed_files =", length(failed_files))
), con = file.path(outdir, "downstream_resume_run_parameters.txt"))

if (length(failed_files) > 0) {
  writeLines(failed_files, con = file.path(outdir, "failed_diffRegulation_files.txt"))
}

cat("\nDownstream analysis finished.\n")
cat("Main directory:", main_dir, "\n")
cat("Output directory:", outdir, "\n")
cat("Failed files:", length(failed_files), "\n")
