S3_BUCKET <- "sod1-app"

S3_BASE_URL <- paste0(
  "https://",
  S3_BUCKET,
  ".s3.ap-southeast-2.amazonaws.com/"
)

get_deg_data <- function(
    dataset,
    level,
    file_path = NULL,
    padj_max = NULL,
    lfc_min = 0,
    direction = "both"
) {
  
  x <- get_deg_data_cloud_rds_cached(
    dataset_id = dataset,
    feature_level = level,
    padj_max = padj_max,
    lfc_min = lfc_min,
    direction = direction
  )
  
  x
}

#local
get_deg_manifest <- function() {
  
  readr::read_csv(
    paste0(S3_BASE_URL, "metadata/deg_manifest.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::mutate(
      path = NA_character_,
      p_str = dplyr::case_when(
        threshold == "padj_0.10" ~ "0.10",
        threshold == "padj_0.05" ~ "0.05",
        threshold == "padj_0.01" ~ "0.01",
        threshold == "padj_0.001" ~ "0.001",
        threshold == "all" ~ NA_character_,
        TRUE ~ NA_character_
      ),
      p = suppressWarnings(as.numeric(p_str))
    ) |>
    dplyr::select(path, dataset, p_str, level, p, threshold)
}

padj_max_to_threshold <- function(padj_max) {
  if (is.null(padj_max) || is.na(padj_max)) {
    return("all")
  }

  dplyr::case_when(
    padj_max == 0.10 ~ "padj_0.10",
    padj_max == 0.05 ~ "padj_0.05",
    padj_max == 0.01 ~ "padj_0.01",
    padj_max == 0.001 ~ "padj_0.001",
    TRUE ~ "all"
  )
}

get_deg_data_cloud_rds <- function(
    dataset_id,
    feature_level,
    padj_max = NULL,
    lfc_min = 0,
    direction = "both"
) {
  threshold <- padj_max_to_threshold(padj_max)

  url <- paste0(
    S3_BASE_URL,
    "deg_results_rds/",
    feature_level, "/",
    dataset_id, "/",
    threshold, ".rds"
  )

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)

  utils::download.file(
    url,
    destfile = tf,
    mode = "wb",
    quiet = TRUE
  )

  x <- readRDS(tf)

  # Apply padj filter in memory when a non-standard threshold (i.e. not one of
  # 0.001, 0.01, 0.05, 0.10) falls back to the "all" prebuilt file.
  if (!is.null(padj_max) && !is.na(padj_max) && threshold == "all" && "padj" %in% names(x)) {
    x <- x[!is.na(x$padj) & x$padj <= padj_max, , drop = FALSE]
  }

  if ("log2FoldChange" %in% names(x)) {
    if (identical(direction, "up")) {
      x <- x[x$log2FoldChange >= lfc_min, , drop = FALSE]
    } else if (identical(direction, "down")) {
      x <- x[x$log2FoldChange <= -lfc_min, , drop = FALSE]
    } else {
      x <- x[abs(x$log2FoldChange) >= lfc_min, , drop = FALSE]
    }
  }

  x
}

deg_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_deg_data_cloud_rds_cached <- memoise::memoise(
  get_deg_data_cloud_rds,
  cache = deg_cache
)


#biomarker server
biomarker_cache <- cachem::cache_mem(max_size = 200 * 1024^2)

get_biomarker_baseline_cloud <- function() {
  url <- paste0(
    S3_BASE_URL,
    "biomarker/baseline_long.rds"
  )

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)

  utils::download.file(url, tf, mode = "wb", quiet = TRUE)
  readRDS(tf)
}

get_biomarker_baseline_cloud_cached <- memoise::memoise(
  get_biomarker_baseline_cloud,
  cache = biomarker_cache
)


#PCA -----------------------------------------------------------
pca_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_pca_manifest_cloud <- function() {
  readr::read_csv(
    paste0(S3_BASE_URL, "metadata/pca_manifest.csv"),
    show_col_types = FALSE
  )
}

get_pca_data_cloud <- function(filename) {

  url <- paste0(
    S3_BASE_URL,
    "pca_input/",
    filename
  )

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)

  utils::download.file(url, tf, mode = "wb", quiet = TRUE)
  readRDS(tf)
}

get_pca_data_cloud_cached <- memoise::memoise(
  get_pca_data_cloud,
  cache = pca_cache
)


# Gene plots
get_tpm_gene_manifest_cloud <- function() {
  readr::read_csv(
    paste0(S3_BASE_URL, "metadata/tpm_gene_manifest.csv"),
    show_col_types = FALSE
  )
}

get_tpm_gene_cloud <- function(filename) {

  url <- paste0(
    S3_BASE_URL,
    "tpm_gene/",
    filename
  )

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  utils::download.file(
    url,
    destfile = tmp,
    mode = "wb",
    quiet = TRUE
  )

  readRDS(tmp)
}

tpm_gene_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_tpm_gene_cloud_cached <- memoise::memoise(
  get_tpm_gene_cloud,
  cache = tpm_gene_cache
)

get_tpm_transcript_manifest_cloud <- function() {
  readr::read_csv(
    paste0(S3_BASE_URL, "metadata/tpm_transcript_manifest.csv"),
    show_col_types = FALSE
  )
}

get_tpm_transcript_cloud <- function(filename) {
  url <- paste0(
    S3_BASE_URL,
    "tpm_transcript/",
    filename
  )

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  utils::download.file(url, tmp, mode = "wb", quiet = TRUE)

  readRDS(tmp)
}

tpm_transcript_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_tpm_transcript_cloud_cached <- memoise::memoise(
  get_tpm_transcript_cloud,
  cache = tpm_transcript_cache
)

get_vsd_manifest_cloud <- function() {
  readr::read_csv(
    paste0(S3_BASE_URL, "metadata/vsd_manifest.csv"),
    show_col_types = FALSE
  )
}


get_vsd_cloud <- function(filename) {
  url <- paste0(
    S3_BASE_URL,
    "vsd/",
    filename
  )

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  utils::download.file(url, tmp, mode = "wb", quiet = TRUE)

  readRDS(tmp)
}

vsd_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_vsd_cloud_cached <- memoise::memoise(
  get_vsd_cloud,
  cache = vsd_cache
)


# GSEA -----------------------------------------------------------

gsea_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_gsea_manifest_cloud <- function() {
  readr::read_csv(
    paste0(S3_BASE_URL, "metadata/gsea_manifest.csv"),
    show_col_types = FALSE
  )
}

get_gsea_data_cloud <- function(filename) {

  url <- paste0(
    S3_BASE_URL,
    "GSEA_results/",
    filename
  )

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)

  utils::download.file(
    url,
    destfile = tf,
    mode = "wb",
    quiet = TRUE
  )

  readRDS(tf)
}

get_gsea_data_cloud_cached <- memoise::memoise(
  get_gsea_data_cloud,
  cache = gsea_cache
)

# GSEA compare / Venn summary ---------------------------------------------

gsea_compare_cache <- cachem::cache_mem(max_size = 100 * 1024^2)

get_gsea_compare_summary_cloud <- function() {
  url <- paste0(
    S3_BASE_URL,
    "GSEA_results/summary/GO_GSEA_direction_summary.csv"
  )

  readr::read_csv(url, show_col_types = FALSE)
}

get_gsea_compare_summary_cloud_cached <- memoise::memoise(
  get_gsea_compare_summary_cloud,
  cache = gsea_compare_cache
)


# GSEA explorer ------------------------------------------------------------

gsea_explore_cache <- cachem::cache_mem(max_size = 200 * 1024^2)

get_gsea_explore_data_cloud <- function() {
  
  url <- paste0(
    S3_BASE_URL,
    "GSEA_results/summary/GSEA_combined_all_datasets.csv"
  )
  
  readr::read_csv(
    url,
    show_col_types = FALSE
  )
}

get_gsea_explore_data_cloud_cached <- memoise::memoise(
  get_gsea_explore_data_cloud,
  cache = gsea_explore_cache
)
