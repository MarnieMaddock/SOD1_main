#' Application Server Logic
#'
#' Defines the server-side logic for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @importFrom ggplot2 ggplot aes geom_point
#' @noRd
app_server <- function(input, output, session) {
  # --- make `pkg` robust for both project + installed package modes ----
  pkg <- tryCatch(utils::packageName(), error = function(e) "")
  if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "SOD1_main"
  pkg <- pkg[[1L]]  # ensure length 1

  pkg_www <- system.file("www", package = pkg, mustWork = FALSE)
  if (nzchar(pkg_www) && dir.exists(pkg_www)) {
    shiny::addResourcePath("pkgwww", pkg_www)  # /pkgwww -> <package>/inst/www
  }

  # Also serve the project copy of inst/www at /projwww (for dev / source tree runs)
  proj_www <- file.path("inst", "www")
  if (dir.exists(proj_www)) {
    shiny::addResourcePath("projwww", proj_www)  # /projwww -> <project>/inst/www
  }

  # --- ensure pretty_map exists (fallback is harmless) ---
  if (!exists("pretty_map", inherits = TRUE)) {
    pretty_map <- stats::setNames(character(0), character(0))
  }

  # --- Call modules (pass the same `pkg`) -------------------------------
  # tpmHeatmapServer("tpm_hm", pkg = pkg, data_mode = data_mode)
  volcanoServer("volc", pkg = pkg)
  GSEAServer("gsea",  pkg = pkg)
  gseaCompareServer("gsea_compare",  pkg = pkg)
  genePlotsServer("gene_plots", pkg = pkg)
  degTablesServer("deg_tables", pkg = pkg)
  degVennServer("deg_venn", pkg = pkg)
  queryGeneAcrossDatasetsServer("gene_query")
  gseaExploreServer("gsea_explore", pkg = pkg)
  # pcaServer("pca", data_mode = data_mode)


}
