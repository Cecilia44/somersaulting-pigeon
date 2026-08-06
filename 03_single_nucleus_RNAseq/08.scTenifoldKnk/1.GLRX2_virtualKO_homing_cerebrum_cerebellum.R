options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(scTenifoldKnk)
  library(future)
  library(furrr)
})

## =========================================
## 1. 用户参数
## =========================================
brain_rds <- "/project-whj/sis/0.temp/2.tumbler/10.SC/3.Seurat+scrublet/3.seurat_DN/fi/cellchat/DN_celltype_res4.rds"
cereb_rds <- "/project-whj/sis/0.temp/2.tumbler/10.SC/3.Seurat+scrublet/3.seurat_XN/fi_2/scRNA_marker_res5_remove50_markerAnnotation.rds"

outdir <- "GLRX2_virtualKO_homing_parallel"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

group_col <- "group"
group_use <- "homing"

brain_celltype_col <- "cell_type"
cereb_celltype_col <- "cell.cluster"

target_candidates <- c("GLRX2", "Glrx2", "glrx2")

exclude_brain_ct <- c(
  "endothelial",
  "vasculature-associated cells",
  "red blood cells"
)

exclude_cereb_ct <- c(
  "red blood cells",
  "endothelial cells",
  "fibroblast",
  "choroid"
)

## -------------------------
## 并行参数
## -------------------------
outer_workers <- 8
inner_ncores  <- 8

## 分析参数
min_cells <- 50
min_pct_expr <- 0.05
max_genes_for_knk <- 2000
padj_cutoff <- 0.05

## 通路库物种
pathway_species <- "Homo sapiens"

## future 全局对象上限
options(future.globals.maxSize = 50 * 1024^3)

## =========================================
## 2. 工具函数
## =========================================
find_target_gene <- function(seu, candidates = c("GLRX2", "Glrx2", "glrx2")) {
  genes <- rownames(seu)

  hit <- candidates[candidates %in% genes]
  if (length(hit) > 0) return(hit[1])

  fuzzy <- genes[grepl("GLRX2|Glrx2|glrx2", genes)]
  if (length(fuzzy) > 0) return(fuzzy[1])

  return(NA_character_)
}

report_glrx2_candidates <- function(seu) {
  genes <- rownames(seu)
  cands1 <- genes[grepl("GLRX|Glrx|glrx", genes)]
  cands2 <- genes[grepl("glutaredoxin|thioredoxin|redox", genes, ignore.case = TRUE)]
  unique(c(cands1, cands2))
}

sanitize_filename <- function(x) {
  x <- gsub("[/ ]", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  x
}

prepare_matrix_for_knk <- function(seu_sub,
                                   target_gene,
                                   min_pct_expr = 0.05,
                                   max_genes_for_knk = 2000) {
  mat <- GetAssayData(seu_sub, assay = "RNA", slot = "counts")

  keep <- Matrix::rowMeans(mat > 0) >= min_pct_expr
  if (target_gene %in% rownames(mat)) {
    keep[rownames(mat) == target_gene] <- TRUE
  }

  mat <- mat[keep, , drop = FALSE]
  mat <- mat[Matrix::rowSums(mat) > 0, , drop = FALSE]

  if (nrow(mat) > max_genes_for_knk) {
    gene_var <- apply(as.matrix(mat), 1, var)
    ord <- order(gene_var, decreasing = TRUE)
    top_genes <- rownames(mat)[ord][seq_len(max_genes_for_knk)]
    top_genes <- unique(c(target_gene, top_genes))
    mat <- mat[intersect(top_genes, rownames(mat)), , drop = FALSE]
  }

  mat
}

run_knk_one_celltype <- function(seu,
                                 tissue_name,
                                 ct,
                                 target_gene,
                                 celltype_col,
                                 min_cells = 50,
                                 min_pct_expr = 0.05,
                                 max_genes_for_knk = 2000,
                                 padj_cutoff = 0.05,
                                 inner_ncores = 8) {

  message("Running: ", tissue_name, " | ", ct)

  cells_use <- rownames(seu@meta.data)[seu@meta.data[[celltype_col]] == ct]
  sub <- subset(seu, cells = cells_use)

  n_cells <- ncol(sub)
  if (n_cells < min_cells) {
    message("  Skip: too few cells (", n_cells, ")")
    return(NULL)
  }

  mat <- prepare_matrix_for_knk(
    seu_sub = sub,
    target_gene = target_gene,
    min_pct_expr = min_pct_expr,
    max_genes_for_knk = max_genes_for_knk
  )

  if (!(target_gene %in% rownames(mat))) {
    message("  Skip: target gene not found after filtering")
    return(NULL)
  }

  res <- tryCatch({
    scTenifoldKnk(
      countMatrix = as.matrix(mat),
      gKO = target_gene,
      qc = TRUE,
      nCores = inner_ncores
    )
  }, error = function(e) {
    message("  ERROR in ", tissue_name, " | ", ct, ": ", e$message)
    return(NULL)
  })

  if (is.null(res) || is.null(res$diffRegulation)) return(NULL)

  dr <- res$diffRegulation %>%
    as.data.frame()

  if (!"gene" %in% colnames(dr)) {
    dr$gene <- rownames(dr)
  }
  if (!"distance" %in% colnames(dr)) dr$distance <- NA_real_
  if (!"p.value" %in% colnames(dr)) dr$p.value <- NA_real_
  if (!"p.adj" %in% colnames(dr)) dr$p.adj <- NA_real_
  if (!"FC" %in% colnames(dr)) dr$FC <- NA_real_
  if (!"Z" %in% colnames(dr)) dr$Z <- NA_real_

  dr <- dr %>%
    mutate(
      tissue = tissue_name,
      cell_type = ct,
      significant = p.adj < padj_cutoff
    ) %>%
    arrange(p.adj, desc(distance))

  summary_df <- tibble(
    tissue = tissue_name,
    cell_type = ct,
    n_cells = n_cells,
    n_genes_input = nrow(mat),
    target_gene = target_gene,
    n_sig = sum(dr$significant, na.rm = TRUE),
    median_distance = median(dr$distance, na.rm = TRUE),
    max_distance = max(dr$distance, na.rm = TRUE),
    mean_distance_sig = ifelse(sum(dr$significant, na.rm = TRUE) > 0,
                               mean(dr$distance[dr$significant], na.rm = TRUE),
                               NA_real_)
  )

  list(
    result = res,
    diff = dr,
    summary = summary_df
  )
}

make_damage_score <- function(summary_tbl) {
  zsafe <- function(x) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    as.numeric(scale(x))
  }

  summary_tbl %>%
    mutate(
      z_n_sig = zsafe(n_sig),
      z_max_distance = zsafe(max_distance),
      z_median_distance = zsafe(median_distance),
      damage_score = z_n_sig + z_max_distance + z_median_distance
    ) %>%
    arrange(desc(damage_score))
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

plot_damage_bar <- function(damage_tbl, file) {
  p <- ggplot(
    damage_tbl,
    aes(x = reorder(paste(tissue, cell_type, sep = " | "), damage_score),
        y = damage_score)
  ) +
    geom_col() +
    coord_flip() +
    labs(
      x = NULL,
      y = "Damage score",
      title = "Predicted vulnerability to GLRX2 virtual knockout (homing only)"
    ) +
    theme_bw(base_size = 12)

  ggsave(file, p, width = 10, height = 8)
}

plot_top_ct <- function(damage_tbl, file, top_n = 20) {
  df <- damage_tbl %>% slice_head(n = top_n)

  p <- ggplot(
    df,
    aes(x = reorder(paste(tissue, cell_type, sep = " | "), damage_score),
        y = damage_score)
  ) +
    geom_col() +
    coord_flip() +
    labs(
      x = NULL,
      y = "Damage score",
      title = paste0("Top ", top_n, " vulnerable cell types")
    ) +
    theme_bw(base_size = 12)

  ggsave(file, p, width = 9, height = 7)
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
    "APOPTOSIS",
    "HYPOXIA",
    "PI3K_AKT_MTOR_SIGNALING",
    "MTORC1_SIGNALING",
    "INFLAMMATORY_RESPONSE",
    "INTERFERON_ALPHA_RESPONSE",
    "INTERFERON_GAMMA_RESPONSE",
    "UNFOLDED_PROTEIN_RESPONSE",
    "MITOCHON",
    "RESPIRATORY",
    "ATP",
    "SYNAPT",
    "CALCIUM",
    "ION"
  )

  fgsea_all %>%
    filter(grepl(paste(key_patterns, collapse = "|"), pathway, ignore.case = TRUE)) %>%
    arrange(database, padj, desc(abs(NES)))
}

## =========================================
## 3. 读入对象
## =========================================
brain <- readRDS(brain_rds)
cereb <- readRDS(cereb_rds)

DefaultAssay(brain) <- "RNA"
DefaultAssay(cereb) <- "RNA"

stopifnot(group_col %in% colnames(brain@meta.data))
stopifnot(group_col %in% colnames(cereb@meta.data))
stopifnot(brain_celltype_col %in% colnames(brain@meta.data))
stopifnot(cereb_celltype_col %in% colnames(cereb@meta.data))

## =========================================
## 4. 筛选 homing
## =========================================
brain_homing <- subset(brain, subset = group == group_use)
cereb_homing <- subset(cereb, subset = group == group_use)

## =========================================
## 5. 过滤不需要的细胞类型
## =========================================
brain_keep_cells <- rownames(brain_homing@meta.data)[
  !(brain_homing@meta.data[[brain_celltype_col]] %in% exclude_brain_ct)
]
brain_homing <- subset(brain_homing, cells = brain_keep_cells)

cereb_keep_cells <- rownames(cereb_homing@meta.data)[
  !(cereb_homing@meta.data[[cereb_celltype_col]] %in% exclude_cereb_ct)
]
cereb_homing <- subset(cereb_homing, cells = cereb_keep_cells)

write.csv(
  as.data.frame(table(brain_homing@meta.data[[brain_celltype_col]])),
  file.path(outdir, "brain_homing_celltype_counts_after_filter.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(cereb_homing@meta.data[[cereb_celltype_col]])),
  file.path(outdir, "cerebellum_homing_celltype_counts_after_filter.csv"),
  row.names = FALSE
)

## =========================================
## 6. 检查 GLRX2
## =========================================
brain_target <- find_target_gene(brain_homing, target_candidates)
cereb_target <- find_target_gene(cereb_homing, target_candidates)

writeLines(c(
  paste("Brain target gene:", brain_target),
  paste("Cerebellum target gene:", cereb_target)
), con = file.path(outdir, "target_gene_detected.txt"))

if (is.na(brain_target)) {
  write.table(
    report_glrx2_candidates(brain_homing),
    file = file.path(outdir, "brain_GLRX2_candidate_genes.txt"),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
}

if (is.na(cereb_target)) {
  write.table(
    report_glrx2_candidates(cereb_homing),
    file = file.path(outdir, "cerebellum_GLRX2_candidate_genes.txt"),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
}

if (is.na(brain_target) | is.na(cereb_target)) {
  stop("GLRX2 not found in one or both filtered objects. Please inspect candidate files.")
}

## =========================================
## 7. 准备通路库
## =========================================
hallmark_df <- msigdbr(species = pathway_species, collection = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

gobp_df <- msigdbr(species = pathway_species, collection = "C5", subcollection = "GO:BP")
gobp_list <- split(gobp_df$gene_symbol, gobp_df$gs_name)

reactome_df <- msigdbr(species = pathway_species, collection = "C2", subcollection = "CP:REACTOME")
reactome_list <- split(reactome_df$gene_symbol, reactome_df$gs_name)

## =========================================
## 8. 设置并行
## =========================================
plan(multisession, workers = outer_workers)

## =========================================
## 9. 分别按细胞类型并行做 virtual KO
## =========================================
brain_cts <- sort(unique(brain_homing@meta.data[[brain_celltype_col]]))
cereb_cts <- sort(unique(cereb_homing@meta.data[[cereb_celltype_col]]))

brain_results <- future_map(
  setNames(brain_cts, brain_cts),
  ~ run_knk_one_celltype(
    seu = brain_homing,
    tissue_name = "brain_homing",
    ct = .x,
    target_gene = brain_target,
    celltype_col = brain_celltype_col,
    min_cells = min_cells,
    min_pct_expr = min_pct_expr,
    max_genes_for_knk = max_genes_for_knk,
    padj_cutoff = padj_cutoff,
    inner_ncores = inner_ncores
  ),
  .options = furrr_options(seed = TRUE)
)

cereb_results <- future_map(
  setNames(cereb_cts, cereb_cts),
  ~ run_knk_one_celltype(
    seu = cereb_homing,
    tissue_name = "cerebellum_homing",
    ct = .x,
    target_gene = cereb_target,
    celltype_col = cereb_celltype_col,
    min_cells = min_cells,
    min_pct_expr = min_pct_expr,
    max_genes_for_knk = max_genes_for_knk,
    padj_cutoff = padj_cutoff,
    inner_ncores = inner_ncores
  ),
  .options = furrr_options(seed = TRUE)
)

all_results <- c(brain_results, cereb_results)
all_results <- all_results[!sapply(all_results, is.null)]

if (length(all_results) == 0) {
  stop("No valid cell types remained after filtering and virtual KO.")
}

plan(sequential)

## =========================================
## 10. 汇总虚拟敲除结果
## =========================================
all_summary <- bind_rows(lapply(all_results, `[[`, "summary"))
all_diff <- bind_rows(lapply(all_results, `[[`, "diff"))

write.csv(
  all_summary,
  file.path(outdir, "GLRX2_virtualKO_homing_celltype_summary_raw.csv"),
  row.names = FALSE
)

write.csv(
  all_diff,
  file.path(outdir, "GLRX2_virtualKO_homing_all_diffRegulation.csv"),
  row.names = FALSE
)

for (nm in names(all_results)) {
  diff_df <- all_results[[nm]]$diff
  tissue_name <- unique(diff_df$tissue)
  ct_name <- unique(diff_df$cell_type)

  fn <- paste0(
    sanitize_filename(tissue_name), "__",
    sanitize_filename(ct_name),
    "__diffRegulation.csv"
  )

  write.csv(diff_df, file.path(outdir, fn), row.names = FALSE)
}

damage_tbl <- make_damage_score(all_summary)

write.csv(
  damage_tbl,
  file.path(outdir, "GLRX2_virtualKO_homing_damage_ranking.csv"),
  row.names = FALSE
)

plot_damage_bar(
  damage_tbl,
  file.path(outdir, "GLRX2_virtualKO_homing_damage_ranking_barplot.pdf")
)

plot_top_ct(
  damage_tbl,
  file.path(outdir, "GLRX2_virtualKO_homing_top20_damage_barplot.pdf"),
  top_n = 20
)

top_gene_tbl <- all_diff %>%
  group_by(tissue, cell_type) %>%
  arrange(p.adj, desc(distance), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  top_gene_tbl,
  file.path(outdir, "GLRX2_virtualKO_homing_top20_perturbed_genes_by_celltype.csv"),
  row.names = FALSE
)

## =========================================
## 11. 对每个细胞类型做通路富集
## =========================================
all_fgsea <- list()

for (nm in names(all_results)) {
  diff_df <- all_results[[nm]]$diff %>%
    filter(!is.na(gene), !is.na(distance))

  if (nrow(diff_df) < 50) next

  tissue_name <- unique(diff_df$tissue)[1]
  ct_name <- unique(diff_df$cell_type)[1]
  base_stub <- paste0(
    sanitize_filename(tissue_name), "__",
    sanitize_filename(ct_name)
  )

  fg_h <- tryCatch({
    safe_fgsea(diff_df, hallmark_list, gene_col = "gene", rank_col = "distance")
  }, error = function(e) NULL)

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

  fg_go <- tryCatch({
    safe_fgsea(diff_df, gobp_list, gene_col = "gene", rank_col = "distance")
  }, error = function(e) NULL)

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

  fg_re <- tryCatch({
    safe_fgsea(diff_df, reactome_list, gene_col = "gene", rank_col = "distance")
  }, error = function(e) NULL)

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
## 12. 记录运行配置
## =========================================
writeLines(c(
  paste("brain_rds =", brain_rds),
  paste("cereb_rds =", cereb_rds),
  paste("group_use =", group_use),
  paste("brain_celltype_col =", brain_celltype_col),
  paste("cereb_celltype_col =", cereb_celltype_col),
  paste("outer_workers =", outer_workers),
  paste("inner_ncores =", inner_ncores),
  paste("min_cells =", min_cells),
  paste("min_pct_expr =", min_pct_expr),
  paste("max_genes_for_knk =", max_genes_for_knk),
  paste("padj_cutoff =", padj_cutoff),
  paste("pathway_species =", pathway_species)
), con = file.path(outdir, "run_parameters.txt"))

cat("\nAnalysis finished.\n")
cat("Output directory:", outdir, "\n")
