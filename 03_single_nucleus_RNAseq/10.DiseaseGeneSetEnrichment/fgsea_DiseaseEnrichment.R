# ============================================================
# fgsea for disease gene sets across cell-type DEGs
# with cell count summary by region / sample / cell type
#
# Major updates:
#   1. Excluded cell types are removed BEFORE FDR correction.
#   2. FDR correction is performed separately within each brain region.
#   3. Dotplots are generated separately for cerebellum and cerebrum.
#   4. Heatmaps are generated separately for cerebellum and cerebrum.
#   5. Heatmaps show NES as color and FDR significance/value as text.
#
# FDR definitions after cell-type exclusion:
#   global_FDR:
#     Within each brain region, adjust all retained cell type x disease gene set p-values.
#
#   FDR_by_celltype:
#     Within each brain region and each retained cell type, adjust p-values
#     across multiple disease gene sets.
#
#   FDR_by_disease:
#     Within each brain region and each disease gene set, adjust p-values
#     across retained cell types.
#
# Dotplot visualization:
#   Point fill color: NES
#   Point size: -log10(FDR)
#   FDR < 0.05: black border
#   0.05 <= FDR < 0.10: dark grey border
#   FDR >= 0.10: light grey border
#
# Heatmap visualization:
#   Tile color: NES
#   ***: FDR < 0.001
#   ** : FDR < 0.01
#   *  : FDR < 0.05
#   FDR >= 0.05: display FDR value with two decimals
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(fgsea)
  library(Seurat)
  library(tidyr)
  library(scales)
})

has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
if (!has_openxlsx) {
  message("Package 'openxlsx' is not installed. Cell-count tables will be saved as CSV instead of XLSX.")
}

# =========================
# 0. Path settings
# =========================

disease_dir <- "./"

deg_dirs <- c(
  cerebellum = "/project-whj/sis/0.temp/2.tumbler/10.SC/3.Seurat+scrublet/3.seurat_XN/fi_2/DEGs",
  cerebrum   = "/project-whj/sis/0.temp/2.tumbler/10.SC/3.Seurat+scrublet/3.seurat_DN/fi/DEG"
)

# Seurat objects for cell count summary
# Modify these paths if file names are different.
seurat_rds_files <- c(
  cerebellum = "/project-whj/sis/0.temp/2.tumbler/10.SC/3.Seurat+scrublet/3.seurat_XN/fi_2/scRNA_marker_res5_remove50_markerAnnotation.rds",
  cerebrum   = "/project-whj/sis/0.temp/2.tumbler/10.SC/3.Seurat+scrublet/3.seurat_DN/fi/cellchat/DN_celltype_res4.rds"
)

# Cell-type metadata column for each brain region
celltype_cols <- c(
  cerebellum = "cell.cluster",
  cerebrum   = "cell_type"
)

# Candidate sample / individual columns
sample_col_candidates <- c("individual.id", "orig.ident", "sample", "sample_id")

# Group column
group_col <- "group"

outdir <- "fgsea_disease_by_celltype"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "gene_sets"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "per_celltype"), showWarnings = FALSE, recursive = TRUE)

# =========================
# 1. Parameters
# =========================

minSize <- 5
maxSize <- 5000
scoreType <- "std"

# Plot all three FDR definitions
fdr_cols_to_plot <- c(
  "global_FDR",
  "FDR_by_celltype",
  "FDR_by_disease"
)

plot_fdr_cutoff <- 0.05
trend_fdr_cutoff <- 0.10

# NES display limit for dotplot and heatmap
nes_limit <- 3

# =========================
# 2. Cell types to exclude BEFORE FDR correction
# =========================

exclude_cerebellum <- c(
  "red blood cells",
  "choroid",
  "fibroblast",
  "endothelial cells"
)

exclude_cerebrum <- c(
  "vasculature-associated cells",
  "endothelial",
  "red blood cells"
)

# =========================
# 3. Helper functions
# =========================

find_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

contains_ignore_case <- function(x, pattern) {
  grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
}

get_first_existing_col <- function(meta, candidates) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

make_fdr_label <- function(fdr) {
  ifelse(
    is.na(fdr), "NA",
    ifelse(
      fdr < 0.001, "***",
      ifelse(
        fdr < 0.01, "**",
        ifelse(
          fdr < 0.05, "*",
          sprintf("%.2f", fdr)
        )
      )
    )
  )
}

# =========================
# 4. Count cells by region / sample / cell type
# =========================

count_cells_by_region_sample_celltype <- function(
  rds_file,
  region_name,
  celltype_col,
  group_col = "group"
) {
  message("Counting cells: ", region_name)

  if (!file.exists(rds_file)) {
    warning("Seurat RDS file not found for ", region_name, ": ", rds_file)
    return(NULL)
  }

  obj <- readRDS(rds_file)
  meta <- as.data.frame(obj@meta.data)

  sample_col <- get_first_existing_col(meta, sample_col_candidates)

  if (is.na(sample_col)) {
    stop(
      "Cannot find sample/individual column in ", region_name,
      ". Tried: ", paste(sample_col_candidates, collapse = ", "),
      "\nAvailable columns:\n", paste(colnames(meta), collapse = ", ")
    )
  }

  if (!celltype_col %in% colnames(meta)) {
    stop(
      "Cannot find cell type column '", celltype_col, "' in ", region_name,
      "\nAvailable columns:\n", paste(colnames(meta), collapse = ", ")
    )
  }

  if (!group_col %in% colnames(meta)) {
    warning("Cannot find group column '", group_col, "' in ", region_name, ". group will be set to NA.")
    meta[[group_col]] <- NA_character_
  }

  count_long <- meta %>%
    mutate(
      region = region_name,
      sample = .data[[sample_col]],
      group = .data[[group_col]],
      cell_type = .data[[celltype_col]]
    ) %>%
    filter(!is.na(sample), !is.na(cell_type), sample != "", cell_type != "") %>%
    count(region, group, sample, cell_type, name = "n_cells") %>%
    arrange(region, group, sample, cell_type)

  count_long
}

cell_count_list <- list()

for (region_i in names(seurat_rds_files)) {
  cell_count_list[[region_i]] <- count_cells_by_region_sample_celltype(
    rds_file = seurat_rds_files[[region_i]],
    region_name = region_i,
    celltype_col = celltype_cols[[region_i]],
    group_col = group_col
  )
}

cell_count_long <- bind_rows(cell_count_list)

if (!is.null(cell_count_long) && nrow(cell_count_long) > 0) {

  cell_count_wide <- cell_count_long %>%
    pivot_wider(
      id_cols = c(region, group, sample),
      names_from = cell_type,
      values_from = n_cells,
      values_fill = 0
    ) %>%
    arrange(region, group, sample)

  cell_count_summary_by_group <- cell_count_long %>%
    group_by(region, group, cell_type) %>%
    summarise(
      total_cells = sum(n_cells),
      n_samples_with_cells = sum(n_cells > 0),
      mean_cells_per_sample = mean(n_cells),
      median_cells_per_sample = median(n_cells),
      min_cells_per_sample = min(n_cells),
      max_cells_per_sample = max(n_cells),
      .groups = "drop"
    ) %>%
    arrange(region, cell_type, group)

  celltype_qc <- cell_count_long %>%
    group_by(region, group, cell_type) %>%
    summarise(
      total_cells = sum(n_cells),
      n_samples_with_cells = sum(n_cells > 0),
      .groups = "drop"
    ) %>%
    pivot_wider(
      id_cols = c(region, cell_type),
      names_from = group,
      values_from = c(total_cells, n_samples_with_cells),
      values_fill = 0
    )

  total_homing_col <- grep("^total_cells_.*homing$", colnames(celltype_qc), value = TRUE)
  total_tumbler_col <- grep("^total_cells_.*tumbler$", colnames(celltype_qc), value = TRUE)
  sample_homing_col <- grep("^n_samples_with_cells_.*homing$", colnames(celltype_qc), value = TRUE)
  sample_tumbler_col <- grep("^n_samples_with_cells_.*tumbler$", colnames(celltype_qc), value = TRUE)

  if (
    length(total_homing_col) == 1 &&
    length(total_tumbler_col) == 1 &&
    length(sample_homing_col) == 1 &&
    length(sample_tumbler_col) == 1
  ) {
    celltype_qc <- celltype_qc %>%
      mutate(
        keep_min50_each_group_sample2 =
          .data[[total_homing_col]] >= 50 &
          .data[[total_tumbler_col]] >= 50 &
          .data[[sample_homing_col]] >= 2 &
          .data[[sample_tumbler_col]] >= 2,

        keep_min100_each_group_sample3 =
          .data[[total_homing_col]] >= 100 &
          .data[[total_tumbler_col]] >= 100 &
          .data[[sample_homing_col]] >= 3 &
          .data[[sample_tumbler_col]] >= 3
      )
  } else {
    warning("Could not identify homing/tumbler columns in celltype_qc. QC flags were not added.")
  }

  if (has_openxlsx) {
    wb <- openxlsx::createWorkbook()

    openxlsx::addWorksheet(wb, "long_counts")
    openxlsx::writeData(wb, "long_counts", cell_count_long)

    openxlsx::addWorksheet(wb, "wide_counts")
    openxlsx::writeData(wb, "wide_counts", cell_count_wide)

    openxlsx::addWorksheet(wb, "summary_by_group")
    openxlsx::writeData(wb, "summary_by_group", cell_count_summary_by_group)

    openxlsx::addWorksheet(wb, "celltype_QC")
    openxlsx::writeData(wb, "celltype_QC", celltype_qc)

    openxlsx::saveWorkbook(
      wb,
      file.path(outdir, "tables", "cell_counts_by_region_sample_celltype.xlsx"),
      overwrite = TRUE
    )
  } else {
    write.csv(
      cell_count_long,
      file.path(outdir, "tables", "cell_counts_long.csv"),
      row.names = FALSE
    )

    write.csv(
      cell_count_wide,
      file.path(outdir, "tables", "cell_counts_wide.csv"),
      row.names = FALSE
    )

    write.csv(
      cell_count_summary_by_group,
      file.path(outdir, "tables", "cell_counts_summary_by_group.csv"),
      row.names = FALSE
    )

    write.csv(
      celltype_qc,
      file.path(outdir, "tables", "celltype_QC_for_fgsea.csv"),
      row.names = FALSE
    )
  }
}

# =========================
# 5. Read PanelApp TSV
# =========================

read_panelapp_green <- function(file) {
  message("Reading PanelApp: ", basename(file))

  dt <- read.delim(
    file = file,
    sep = "\t",
    header = TRUE,
    quote = "\"",
    fill = TRUE,
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = TRUE
  )

  colnames(dt) <- make.names(colnames(dt))

  entity_col <- find_first_existing_col(dt, c("Entity.type", "Entity.Type", "Entity_type"))
  gene_col   <- find_first_existing_col(dt, c("Gene.Symbol", "Gene.symbol", "Gene"))
  gel_col    <- find_first_existing_col(dt, c("GEL_Status", "GEL.status"))
  source_col <- find_first_existing_col(dt, c("Sources...separated.", "Sources"))

  if (is.na(entity_col)) {
    stop(
      "Cannot find entity column in ", basename(file),
      "\nAvailable columns:\n", paste(colnames(dt), collapse = ", ")
    )
  }

  if (is.na(gene_col)) {
    stop(
      "Cannot find gene symbol column in ", basename(file),
      "\nAvailable columns:\n", paste(colnames(dt), collapse = ", ")
    )
  }

  dt2 <- dt[
    dt[[entity_col]] == "gene" &
      !is.na(dt[[gene_col]]) &
      dt[[gene_col]] != "",
    ,
    drop = FALSE
  ]

  if (!is.na(gel_col)) {
    suppressWarnings(dt2[[gel_col]] <- as.numeric(dt2[[gel_col]]))
    dt2 <- dt2[dt2[[gel_col]] == 3, , drop = FALSE]
  } else if (!is.na(source_col)) {
    dt2 <- dt2[contains_ignore_case(dt2[[source_col]], "green"), , drop = FALSE]
  } else {
    warning("No GEL_Status/Sources column found in ", basename(file), "; returning all gene entities.")
  }

  genes <- sort(unique(trimws(dt2[[gene_col]])))
  genes <- genes[!is.na(genes) & genes != ""]

  genes
}

# =========================
# 6. Read ClinGen
# =========================

read_clingen <- function(file) {
  message("Reading ClinGen: ", basename(file))

  dt <- fread(
    file = file,
    skip = 4,
    data.table = FALSE
  )

  colnames(dt) <- trimws(colnames(dt))

  required_cols <- c("GENE SYMBOL", "DISEASE LABEL", "CLASSIFICATION")
  miss <- setdiff(required_cols, colnames(dt))

  if (length(miss) > 0) {
    stop(
      "ClinGen file missing columns: ", paste(miss, collapse = ", "),
      "\nAvailable columns:\n", paste(colnames(dt), collapse = ", ")
    )
  }

  dt <- dt %>%
    filter(!is.na(`GENE SYMBOL`), `GENE SYMBOL` != "", `GENE SYMBOL` != "++++++++++++") %>%
    mutate(
      disease_lower = tolower(`DISEASE LABEL`),
      class_upper   = toupper(CLASSIFICATION)
    )

  dt
}

# =========================
# 7. Build disease gene sets
# =========================

f_epi <- file.path(disease_dir, "Early onset or syndromic epilepsy.tsv")

f_ataxia_adult <- file.path(disease_dir, "Hereditary ataxia with onset in adulthood.tsv")
f_ataxia_child <- file.path(disease_dir, "Hereditary ataxia and cerebellar anomalies - childhood onset.tsv")

f_hsp_adult <- file.path(disease_dir, "Adult onset hereditary spastic paraplegia.tsv")
f_hsp_child <- file.path(disease_dir, "Childhood onset hereditary spastic paraplegia.tsv")

f_dystonia_adult <- file.path(disease_dir, "Adult onset dystonia, chorea or related movement disorder.tsv")
f_dystonia_child <- file.path(disease_dir, "Childhood onset dystonia, chorea or related movement disorder.tsv")

f_clingen <- file.path(disease_dir, "Clingen-Gene-Disease-Summary-2026-04-09(1).csv")
if (!file.exists(f_clingen)) {
  f_clingen <- file.path(disease_dir, "Clingen-Gene-Disease-Summary-2026-04-09.csv")
}

input_files <- c(
  f_epi,
  f_ataxia_adult,
  f_ataxia_child,
  f_hsp_adult,
  f_hsp_child,
  f_dystonia_adult,
  f_dystonia_child,
  f_clingen
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("These input files were not found:\n", paste(missing_files, collapse = "\n"))
}

epi_genes <- read_panelapp_green(f_epi)

ataxia_adult_genes <- read_panelapp_green(f_ataxia_adult)
ataxia_child_genes <- read_panelapp_green(f_ataxia_child)
ataxia_genes <- sort(unique(c(ataxia_adult_genes, ataxia_child_genes)))

hsp_adult_genes <- read_panelapp_green(f_hsp_adult)
hsp_child_genes <- read_panelapp_green(f_hsp_child)
hsp_genes <- sort(unique(c(hsp_adult_genes, hsp_child_genes)))

dystonia_adult_genes <- read_panelapp_green(f_dystonia_adult)
dystonia_child_genes <- read_panelapp_green(f_dystonia_child)
dystonia_genes <- sort(unique(c(dystonia_adult_genes, dystonia_child_genes)))

clingen <- read_clingen(f_clingen)

pd_genes <- clingen %>%
  filter(grepl("parkinson", disease_lower, ignore.case = TRUE)) %>%
  filter(class_upper %in% c("DEFINITIVE")) %>%
  pull(`GENE SYMBOL`) %>%
  unique() %>%
  sort()

mito_genes <- clingen %>%
  filter(grepl("mitochond", disease_lower, ignore.case = TRUE)) %>%
  filter(class_upper %in% c("DEFINITIVE", "STRONG")) %>%
  pull(`GENE SYMBOL`) %>%
  unique() %>%
  sort()

disease_sets <- list(
  epilepsy              = epi_genes,
  ataxia                = ataxia_genes,
  Parkinsons_disease    = pd_genes,
  HSP                   = hsp_genes,
  dystonia              = dystonia_genes,
  mitochondrial_disease = mito_genes
)

disease_sets <- disease_sets[sapply(disease_sets, length) > 0]

for (nm in names(disease_sets)) {
  writeLines(disease_sets[[nm]], file.path(outdir, "gene_sets", paste0(nm, ".txt")))
}

gene_set_sizes <- data.frame(
  disease = names(disease_sets),
  n_genes = sapply(disease_sets, length),
  stringsAsFactors = FALSE
)

write.csv(
  gene_set_sizes,
  file.path(outdir, "tables", "disease_gene_set_sizes.csv"),
  row.names = FALSE
)

print(gene_set_sizes)

# =========================
# 8. Read DEG files
# =========================

get_deg_files <- function(dir_path, region_name) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  data.frame(file = files, region = region_name, stringsAsFactors = FALSE)
}

deg_file_info <- rbind(
  get_deg_files(deg_dirs["cerebellum"], "cerebellum"),
  get_deg_files(deg_dirs["cerebrum"],   "cerebrum")
)

deg_file_info <- deg_file_info[
  grepl("degs\\.csv$", basename(deg_file_info$file)),
  ,
  drop = FALSE
]

if (nrow(deg_file_info) == 0) {
  stop("No DEG files matching '*degs.csv' were found.")
}

extract_celltype_from_filename <- function(x) {
  x <- basename(x)
  x <- sub("_tumbler_vs_homing_degs\\.csv$", "", x)
  x
}

read_deg_file <- function(file, region) {
  dt <- fread(file = file, data.table = FALSE)
  colnames(dt)[1] <- "gene"

  required_cols <- c("gene", "stat", "log2FoldChange", "padj")
  miss <- setdiff(required_cols, colnames(dt))

  if (length(miss) > 0) {
    stop(
      "DEG file missing columns in ", basename(file), ": ",
      paste(miss, collapse = ", ")
    )
  }

  dt <- dt %>%
    mutate(
      gene = trimws(gene),
      cell_type = extract_celltype_from_filename(file),
      region = region
    ) %>%
    filter(!is.na(gene), gene != "", !is.na(stat))

  dt
}

deg_list <- Map(read_deg_file, deg_file_info$file, deg_file_info$region)
deg_all <- do.call(rbind, deg_list)

write.csv(
  deg_all,
  file.path(outdir, "tables", "all_deg_merged.csv"),
  row.names = FALSE
)

# =========================
# 9. Prepare preranked stats
# =========================

prepare_ranked_stats <- function(df) {
  df2 <- df %>%
    filter(!is.na(gene), gene != "", !is.na(stat)) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()

  ranks <- df2$stat
  names(ranks) <- df2$gene
  ranks <- sort(ranks, decreasing = TRUE)

  ranks
}

# =========================
# 10. Run fgsea for one cell type
# =========================

run_fgsea_one <- function(
  df_sub,
  disease_sets,
  minSize = 5,
  maxSize = 5000,
  scoreType = "std"
) {
  ranks <- prepare_ranked_stats(df_sub)

  pathways_use <- lapply(disease_sets, function(gs) intersect(gs, names(ranks)))
  pathways_use <- pathways_use[sapply(pathways_use, length) >= minSize]

  if (length(pathways_use) == 0) {
    return(data.frame())
  }

  fg <- suppressWarnings(
    fgseaMultilevel(
      pathways = pathways_use,
      stats = ranks,
      minSize = minSize,
      maxSize = maxSize,
      scoreType = scoreType
    )
  )

  fg <- as.data.frame(fg)
  if (nrow(fg) == 0) return(data.frame())

  if ("leadingEdge" %in% colnames(fg)) {
    fg$leadingEdge_genes <- vapply(
      fg$leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
    fg$leadingEdge <- NULL
  }

  fg$pathway_size_in_universe <- vapply(
    pathways_use[fg$pathway],
    length,
    integer(1)
  )

  fg
}

# =========================
# 11. Run fgsea for all region x cell type
# =========================

group_keys <- deg_all %>%
  distinct(region, cell_type)

fgsea_res_list <- vector("list", nrow(group_keys))

for (i in seq_len(nrow(group_keys))) {
  region_i <- group_keys$region[i]
  cell_i   <- group_keys$cell_type[i]

  message("Running fgsea: ", region_i, " | ", cell_i)

  sub <- deg_all %>%
    filter(region == region_i, cell_type == cell_i)

  fg <- run_fgsea_one(
    df_sub = sub,
    disease_sets = disease_sets,
    minSize = minSize,
    maxSize = maxSize,
    scoreType = scoreType
  )

  if (nrow(fg) == 0) next

  fg$region <- region_i
  fg$cell_type <- cell_i
  fg$n_genes_ranked <- length(unique(sub$gene))

  is_list_col <- vapply(fg, is.list, logical(1))
  if (any(is_list_col)) {
    fg <- fg[, !is_list_col, drop = FALSE]
  }

  out_file <- file.path(
    outdir,
    "per_celltype",
    paste0(region_i, "__", cell_i, "__fgsea.csv")
  )
  write.csv(fg, out_file, row.names = FALSE)

  fgsea_res_list[[i]] <- fg
}

fgsea_res <- do.call(rbind, fgsea_res_list)

if (is.null(fgsea_res) || nrow(fgsea_res) == 0) {
  stop("fgsea returned no results. Check minSize and gene symbol overlap.")
}

is_list_col <- vapply(fgsea_res, is.list, logical(1))
if (any(is_list_col)) {
  fgsea_res <- fgsea_res[, !is_list_col, drop = FALSE]
}

# =========================
# 12. Remove excluded cell types BEFORE FDR correction
# =========================
# Important:
#   These cell types are removed before global_FDR,
#   FDR_by_celltype, and FDR_by_disease are calculated.
#   Therefore, they will not affect multiple-testing correction.

fgsea_res_raw_with_excluded <- fgsea_res

write.csv(
  fgsea_res_raw_with_excluded,
  file.path(outdir, "tables", "fgsea_disease_all_celltypes_BEFORE_exclusion.csv"),
  row.names = FALSE
)

fgsea_res <- fgsea_res %>%
  filter(
    !(region == "cerebellum" & cell_type %in% exclude_cerebellum),
    !(region == "cerebrum" & cell_type %in% exclude_cerebrum)
  )

write.csv(
  fgsea_res,
  file.path(outdir, "tables", "fgsea_disease_all_celltypes_AFTER_exclusion_before_FDR.csv"),
  row.names = FALSE
)

excluded_rows <- fgsea_res_raw_with_excluded %>%
  anti_join(
    fgsea_res %>% select(region, cell_type, pathway),
    by = c("region", "cell_type", "pathway")
  )

write.csv(
  excluded_rows,
  file.path(outdir, "tables", "fgsea_disease_excluded_celltypes_rows.csv"),
  row.names = FALSE
)

message("Cell types excluded before FDR correction:")
message("  cerebellum: ", paste(exclude_cerebellum, collapse = ", "))
message("  cerebrum: ", paste(exclude_cerebrum, collapse = ", "))
message("Rows before exclusion: ", nrow(fgsea_res_raw_with_excluded))
message("Rows after exclusion: ", nrow(fgsea_res))
message("Rows excluded: ", nrow(excluded_rows))

# =========================
# 13. Region-separated FDR correction after exclusion
# =========================

fgsea_res <- fgsea_res %>%
  group_by(region) %>%
  mutate(global_FDR = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  group_by(region, cell_type) %>%
  mutate(FDR_by_celltype = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  group_by(region, pathway) %>%
  mutate(FDR_by_disease = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  arrange(region, FDR_by_celltype, global_FDR, FDR_by_disease, padj, pathway)

write.csv(
  fgsea_res,
  file.path(outdir, "tables", "fgsea_disease_all_celltypes.csv"),
  row.names = FALSE
)

write.csv(
  fgsea_res %>% filter(global_FDR < 0.05),
  file.path(outdir, "tables", "fgsea_disease_significant_globalFDR.csv"),
  row.names = FALSE
)

write.csv(
  fgsea_res %>% filter(FDR_by_celltype < 0.05),
  file.path(outdir, "tables", "fgsea_disease_significant_byCelltypeFDR.csv"),
  row.names = FALSE
)

write.csv(
  fgsea_res %>% filter(FDR_by_disease < 0.05),
  file.path(outdir, "tables", "fgsea_disease_significant_byDiseaseFDR.csv"),
  row.names = FALSE
)

write.csv(
  fgsea_res %>% filter(FDR_by_celltype < 0.10),
  file.path(outdir, "tables", "fgsea_disease_trend_byCelltypeFDR_0.10.csv"),
  row.names = FALSE
)

# Region-specific tables
for (region_i in unique(fgsea_res$region)) {

  sub_region <- fgsea_res %>%
    filter(region == region_i)

  write.csv(
    sub_region,
    file.path(outdir, "tables", paste0("fgsea_disease_all_celltypes_", region_i, ".csv")),
    row.names = FALSE
  )

  write.csv(
    sub_region %>% filter(global_FDR < 0.05),
    file.path(outdir, "tables", paste0("fgsea_disease_significant_globalFDR_", region_i, ".csv")),
    row.names = FALSE
  )

  write.csv(
    sub_region %>% filter(FDR_by_celltype < 0.05),
    file.path(outdir, "tables", paste0("fgsea_disease_significant_byCelltypeFDR_", region_i, ".csv")),
    row.names = FALSE
  )

  write.csv(
    sub_region %>% filter(FDR_by_celltype < 0.10),
    file.path(outdir, "tables", paste0("fgsea_disease_trend_byCelltypeFDR_0.10_", region_i, ".csv")),
    row.names = FALSE
  )

  write.csv(
    sub_region %>% filter(FDR_by_disease < 0.05),
    file.path(outdir, "tables", paste0("fgsea_disease_significant_byDiseaseFDR_", region_i, ".csv")),
    row.names = FALSE
  )
}

# =========================
# 14. Prepare plotting table
# =========================

plot_df <- fgsea_res %>%
  mutate(
    cell_label = paste(region, cell_type, sep = " | "),
    neglog10_globalFDR = -log10(global_FDR + 1e-300),
    neglog10_FDR_by_celltype = -log10(FDR_by_celltype + 1e-300),
    neglog10_FDR_by_disease = -log10(FDR_by_disease + 1e-300),
    direction = ifelse(
      NES > 0,
      "Higher in tumbler",
      ifelse(NES < 0, "Lower in tumbler", "NA")
    )
  ) %>%
  mutate(
    sig_mark_global = ifelse(
      global_FDR < 0.001, "***",
      ifelse(global_FDR < 0.01, "**",
             ifelse(global_FDR < 0.05, "*", ""))
    ),
    sig_mark_by_celltype = ifelse(
      FDR_by_celltype < 0.001, "***",
      ifelse(FDR_by_celltype < 0.01, "**",
             ifelse(FDR_by_celltype < 0.05, "*", ""))
    ),
    sig_mark_by_disease = ifelse(
      FDR_by_disease < 0.001, "***",
      ifelse(FDR_by_disease < 0.01, "**",
             ifelse(FDR_by_disease < 0.05, "*", ""))
    )
  )

write.csv(
  plot_df,
  file.path(outdir, "tables", "fgsea_disease_filtered_celltypes.csv"),
  row.names = FALSE
)

for (region_i in unique(plot_df$region)) {
  write.csv(
    plot_df %>% filter(region == region_i),
    file.path(outdir, "tables", paste0("fgsea_disease_filtered_celltypes_", region_i, ".csv")),
    row.names = FALSE
  )
}

# =========================
# 15. Cerebellum Golgi-specific outputs
# =========================

golgi_res <- fgsea_res %>%
  filter(region == "cerebellum", tolower(cell_type) == "golgi") %>%
  arrange(FDR_by_celltype, pval)

if (nrow(golgi_res) > 0) {
  write.csv(
    golgi_res,
    file.path(outdir, "tables", "cerebellum_golgi_fgsea_disease_results.csv"),
    row.names = FALSE
  )

  write.csv(
    golgi_res %>% filter(pval < 0.05 | FDR_by_celltype < 0.10),
    file.path(outdir, "tables", "cerebellum_golgi_fgsea_disease_trend_results.csv"),
    row.names = FALSE
  )
}

golgi_deg <- deg_all %>%
  filter(region == "cerebellum", tolower(cell_type) == "golgi")

if (nrow(golgi_deg) > 0) {
  golgi_ranks <- prepare_ranked_stats(golgi_deg)

  golgi_overlap <- data.frame(
    disease = names(disease_sets),
    disease_set_size = sapply(disease_sets, length),
    overlap_with_golgi_ranked_genes = sapply(disease_sets, function(gs) {
      length(intersect(gs, names(golgi_ranks)))
    }),
    stringsAsFactors = FALSE
  )

  write.csv(
    golgi_overlap,
    file.path(outdir, "tables", "cerebellum_golgi_disease_gene_overlap.csv"),
    row.names = FALSE
  )
}

# =========================
# 16. Annotated NES heatmap
# =========================

plot_annotated_nes_heatmap <- function(
  df,
  region_i,
  fdr_col,
  outdir,
  file_prefix = "fgsea_NES_annotated_heatmap",
  width = 8.8,
  height = NULL,
  nes_limit = 3
) {

  df_plot <- df %>%
    filter(region == region_i) %>%
    mutate(
      plot_FDR = .data[[fdr_col]],
      NES_plot = pmax(pmin(NES, nes_limit), -nes_limit),
      fdr_label = make_fdr_label(plot_FDR),
      cell_label_clean = cell_type,
      text_color = ifelse(abs(NES_plot) >= 2.2, "white", "black")
    )

  if (nrow(df_plot) == 0) return(NULL)

  cell_levels <- rev(unique(df_plot$cell_label_clean))
  pathway_levels <- unique(df_plot$pathway)

  df_plot <- df_plot %>%
    mutate(
      cell_label_clean = factor(cell_label_clean, levels = cell_levels),
      pathway = factor(pathway, levels = pathway_levels)
    )

  if (is.null(height)) {
    height <- max(5.5, 0.45 * length(cell_levels))
  }

  p_heat <- ggplot(df_plot, aes(x = pathway, y = cell_label_clean, fill = NES_plot)) +
    geom_tile(color = "grey85", linewidth = 0.25) +
    geom_text(
      aes(label = fdr_label, color = text_color),
      size = 3.2,
      fontface = "bold"
    ) +
    scale_color_identity() +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(-nes_limit, nes_limit),
      oob = scales::squish,
      name = "NES"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(
        "Disease gene set fgsea NES in ",
        region_i,
        "; significance by ",
        fdr_col
      ),
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      plot.margin = margin(10, 18, 10, 10)
    )

  ggsave(
    file.path(
      outdir,
      "plots",
      paste0(file_prefix, "_", region_i, "__", fdr_col, ".pdf")
    ),
    p_heat,
    width = width,
    height = height,
    limitsize = FALSE
  )

  ggsave(
    file.path(
      outdir,
      "plots",
      paste0(file_prefix, "_", region_i, "__", fdr_col, ".png")
    ),
    p_heat,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE
  )

  return(p_heat)
}

# =========================
# 17. Plot one FDR set
# =========================

plot_one_fdr_set <- function(
  plot_df_input,
  fdr_col,
  fdr_cutoff = 0.05,
  trend_cutoff = 0.10,
  nes_limit = 3
) {

  message("Plotting figures using: ", fdr_col)

  df <- plot_df_input

  df$plot_FDR <- df[[fdr_col]]
  df$neglog10_plotFDR <- -log10(df$plot_FDR + 1e-300)
  df$sig_plot <- !is.na(df$plot_FDR) & df$plot_FDR < fdr_cutoff

  df <- df %>%
    mutate(
      sig_mark = ifelse(
        plot_FDR < 0.001, "***",
        ifelse(plot_FDR < 0.01, "**",
               ifelse(plot_FDR < 0.05, "*", ""))
      )
    )

  write.csv(
    df,
    file.path(outdir, "tables", paste0("fgsea_disease_plot_table__", fdr_col, ".csv")),
    row.names = FALSE
  )

  for (region_i in unique(df$region)) {
    write.csv(
      df %>% filter(region == region_i),
      file.path(outdir, "tables", paste0("fgsea_disease_plot_table_", region_i, "__", fdr_col, ".csv")),
      row.names = FALSE
    )
  }

  # -------------------------
  # Region-separated dotplots
  # -------------------------

  for (region_i in unique(df$region)) {

    df_region <- df %>%
      filter(region == region_i)

    if (nrow(df_region) == 0) next

    df_region$cell_label <- factor(
      df_region$cell_label,
      levels = rev(unique(df_region$cell_label))
    )

    df_region <- df_region %>%
      mutate(
        NES_plot = pmax(pmin(NES, nes_limit), -nes_limit),
        logFDR_plot = pmax(neglog10_plotFDR, 0),
        sig_group = case_when(
          !is.na(plot_FDR) & plot_FDR < fdr_cutoff ~ paste0(fdr_col, " < ", fdr_cutoff),
          !is.na(plot_FDR) & plot_FDR < trend_cutoff ~ paste0(fdr_cutoff, " <= ", fdr_col, " < ", trend_cutoff),
          TRUE ~ paste0(fdr_col, " >= ", trend_cutoff)
        )
      )

    df_ns <- df_region %>%
      filter(is.na(plot_FDR) | plot_FDR >= trend_cutoff)

    df_trend <- df_region %>%
      filter(!is.na(plot_FDR), plot_FDR >= fdr_cutoff, plot_FDR < trend_cutoff)

    df_sig <- df_region %>%
      filter(!is.na(plot_FDR), plot_FDR < fdr_cutoff)

    p_dot_region <- ggplot() +
      geom_point(
        data = df_ns,
        aes(x = pathway, y = cell_label, size = logFDR_plot, fill = NES_plot),
        shape = 21,
        color = "grey72",
        stroke = 0.25,
        alpha = 0.80
      ) +
      geom_point(
        data = df_trend,
        aes(x = pathway, y = cell_label, size = logFDR_plot, fill = NES_plot),
        shape = 21,
        color = "grey35",
        stroke = 0.55,
        alpha = 0.95
      ) +
      geom_point(
        data = df_sig,
        aes(x = pathway, y = cell_label, size = logFDR_plot, fill = NES_plot),
        shape = 21,
        color = "black",
        stroke = 0.85,
        alpha = 0.98
      ) +
      scale_fill_gradient2(
        low = "#3B4CC0",
        mid = "white",
        high = "#B40426",
        midpoint = 0,
        limits = c(-nes_limit, nes_limit),
        oob = scales::squish,
        name = "NES"
      ) +
      scale_size_continuous(
        range = c(3, 10),
        name = paste0("-log10(", fdr_col, ")")
      ) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(
          "Disease gene set fgsea in ",
          region_i,
          "; significance by ",
          fdr_col
        ),
        x = NULL,
        y = NULL
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        plot.margin = margin(10, 18, 10, 10)
      )

    ggsave(
      file.path(
        outdir,
        "plots",
        paste0("fgsea_disease_dotplot_", region_i, "__", fdr_col, ".pdf")
      ),
      p_dot_region,
      width = 9.5,
      height = max(5.5, 0.42 * length(unique(df_region$cell_label))),
      limitsize = FALSE
    )

    ggsave(
      file.path(
        outdir,
        "plots",
        paste0("fgsea_disease_dotplot_", region_i, "__", fdr_col, ".png")
      ),
      p_dot_region,
      width = 9.5,
      height = max(5.5, 0.42 * length(unique(df_region$cell_label))),
      dpi = 300,
      limitsize = FALSE
    )
  }

  # -------------------------
  # Region-separated annotated NES heatmaps
  # -------------------------

  for (region_i in unique(df$region)) {

    df_region <- df %>%
      filter(region == region_i)

    if (nrow(df_region) == 0) next

    plot_annotated_nes_heatmap(
      df = df_region,
      region_i = region_i,
      fdr_col = fdr_col,
      outdir = outdir,
      file_prefix = "fgsea_NES_annotated_heatmap",
      width = 8.8,
      height = max(5.5, 0.45 * length(unique(df_region$cell_type))),
      nes_limit = nes_limit
    )
  }
}

# =========================
# 18. Generate plots for all FDR definitions
# =========================

for (fdr_i in fdr_cols_to_plot) {
  plot_one_fdr_set(
    plot_df_input = plot_df,
    fdr_col = fdr_i,
    fdr_cutoff = plot_fdr_cutoff,
    trend_cutoff = trend_fdr_cutoff,
    nes_limit = nes_limit
  )
}

# =========================
# 19. Main messages
# =========================

message("All done.")
message("Results saved in: ", normalizePath(outdir))

message("Important filtering strategy:")
message("  Excluded cell types were removed BEFORE FDR correction.")
message("  Therefore, they do not affect global_FDR or FDR_by_disease.")

message("Excluded cell types:")
message("  cerebellum: ", paste(exclude_cerebellum, collapse = ", "))
message("  cerebrum: ", paste(exclude_cerebrum, collapse = ", "))

message("FDR correction strategy:")
message("  global_FDR: adjusted within each brain region after cell-type exclusion")
message("  FDR_by_celltype: adjusted within each brain region and each retained cell type")
message("  FDR_by_disease: adjusted within each brain region and each disease gene set across retained cell types")

message("Main FDR_by_celltype dotplots:")
message("  ", file.path(outdir, "plots", "fgsea_disease_dotplot_cerebellum__FDR_by_celltype.pdf"))
message("  ", file.path(outdir, "plots", "fgsea_disease_dotplot_cerebrum__FDR_by_celltype.pdf"))

message("Main FDR_by_celltype annotated heatmaps:")
message("  ", file.path(outdir, "plots", "fgsea_NES_annotated_heatmap_cerebellum__FDR_by_celltype.pdf"))
message("  ", file.path(outdir, "plots", "fgsea_NES_annotated_heatmap_cerebrum__FDR_by_celltype.pdf"))

message("Diagnostic tables:")
message("  ", file.path(outdir, "tables", "fgsea_disease_all_celltypes_BEFORE_exclusion.csv"))
message("  ", file.path(outdir, "tables", "fgsea_disease_all_celltypes_AFTER_exclusion_before_FDR.csv"))
message("  ", file.path(outdir, "tables", "fgsea_disease_excluded_celltypes_rows.csv"))
