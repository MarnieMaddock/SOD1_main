#' @importFrom shiny HTML
#' @noRd
genePlotsSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Gene plots"),
    checkboxGroupInput(
      ns("gp_datasets"),
      label   = "Datasets",
      choices = character(0)   # filled by server
    ),
    uiOutput(ns("gp_datasets_note")),
    tags$br(),
    textInput(
      ns("gp_gene"),
      "Gene",
      value = "SOD1",   # pre-filled
      placeholder = "Type a gene symbol to search..."
    ),
    tags$br(),
    checkboxInput(ns("gp_logy"), "Log10 Y-axis", value = FALSE),
    tags$hr(),
    strong("Download"),
    br(),
    downloadButton(ns("dl_points"), "Replicates (CSV)"),
    downloadButton(ns("dl_summary"), "Summary (CSV)"),
    downloadButton(ns("dl_plot"), "Plot (SVG)"),
    downloadButton(ns("dl_plot_png"), "Plot (PNG)")
  )
}

#' Gene Plots - main UI
#' @noRd
genePlotsMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    # if you don't have shinycssloaders installed yet, use plotOutput directly
    if (requireNamespace("shinycssloaders", quietly = TRUE)) {
      shinycssloaders::withSpinner(plotOutput(ns("gp_plot"), height = "60vh", width = "90vh"), type = 4,  color = "#005249")
    } else {
      plotOutput(ns("gp_plot"), height = "60vh", width= "90vh")
    },
    br(),
    DT::dataTableOutput(ns("gp_table"))
  )
}

#' Gene Plots server - multi-dataset (checkboxes)
#' @noRd
genePlotsServer <- function(id, pkg = "SOD1main"){

  moduleServer(id, function(input, output, session) {

    pkg <- tryCatch(pkg, error = function(e) "")
    if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) {
      pkg <- "SOD1main"
    }
    pkg <- pkg[[1L]]
    pretty_map <- c(
      "all_data" = "SOD1 knockout"
    )

    `%||%` <- function(a, b) if (is.null(a)) b else a
    pretty_label <- function(id) pretty_map[[id]] %||% id

    if (!exists("theme_Marnie", inherits = TRUE)) {
      tp <- system.file("R", "utils_graphTheme.R", package = pkg, mustWork = FALSE)
      if (!nzchar(tp)) tp <- file.path("R", "utils_graphTheme.R")
      if (file.exists(tp)) source(tp)
    }

    tpm_manifest <- reactiveVal(NULL)
    
    observeEvent(TRUE, {
      
      man <- get_tpm_gene_manifest_cloud()
      
      validate(
        need(
          "filename" %in% names(man),
          "TPM manifest needs a filename column."
        ),
        need(
          "dataset" %in% names(man),
          "TPM manifest needs a dataset column."
        )
      )
      
      tpm_manifest(man)
      
    }, once = TRUE)

    manifest <- reactive({
      req(tpm_manifest())
      tpm_manifest()
    })

    observe({
      m <- manifest()

      ds <- sort(intersect(unique(m$dataset), names(pretty_map)))

      validate(
        need(length(ds) > 0, "No TPM gene datasets found.")
      )

      choices_named <- stats::setNames(ds, unname(pretty_map[ds]))
      pre_sel <- ds[1]

      updateCheckboxGroupInput(
        session,
        "gp_datasets",
        choices = choices_named,
        selected = pre_sel
      )
    })



    output$gp_datasets_note <- renderUI({
      req(input$gp_datasets)

      labs <- vapply(
        input$gp_datasets,
        pretty_label,
        "",
        USE.NAMES = FALSE
      )

      htmltools::HTML(
        sprintf(
          "<small>Selected: <b>%s</b></small>",
          paste(labs, collapse = ", ")
        )
      )
    })

    file_paths <- reactive({
      
      req(input$gp_datasets)
      
      m <- manifest() |>
        dplyr::filter(dataset %in% input$gp_datasets)
      
      validate(
        need(
          nrow(m) >= 1,
          "No TPM files found for the selected datasets."
        )
      )
      
      stats::setNames(
        m$filename,
        m$dataset
      )
    })

    tpms_wide_list <- reactive({
      
      fps <- file_paths()
      
      out <- lapply(names(fps), function(ds) {
        
        x <- get_tpm_gene_cloud_cached(
          fps[[ds]]
        )
        
        if (!is.data.frame(x)) {
          x <- as.data.frame(x)
        }
        
        validate(
          need(
            all(c("gene_id", "gene_name") %in% names(x)),
            sprintf(
              "TPM RDS for '%s' must contain 'gene_id' and 'gene_name'.",
              ds
            )
          )
        )
        
        x
      })
      
      names(out) <- names(fps)
      
      out
    })

    observeEvent(tpms_wide_list(), {
      wl <- tpms_wide_list()

      genes <- sort(unique(unlist(lapply(wl, function(x) x$gene_name))))

      validate(
        need(length(genes) > 0, "No gene names found in the selected TPM files.")
      )

      current <- isolate(input$gp_gene)

      selected <- if (!is.null(current) && nzchar(current) && current %in% genes) {
        current
      } else if ("SOD1" %in% genes) {
        "SOD1"
      } else {
        genes[1]
      }

      updateTextInput(
        session,
        "gp_gene",
        value = selected
      )

    }, ignoreInit = FALSE)

    tpms_long_all <- reactive({
      wl <- tpms_wide_list()
      if (!length(wl)) return(dplyr::tibble())
      parts <- lapply(names(wl), function(ds) {
        x <- wl[[ds]]
        lng <- tidyr::pivot_longer(
          x, cols = -c(gene_id, gene_name),
          names_to = "sample", values_to = "TPM"
        )
        lng$dataset <- ds
        lng
      })
      dplyr::bind_rows(parts)
    })


    parse_conditions <- function(df_ds) {
      
      repnum <- suppressWarnings(
        as.integer(
          sub(".*_Rep([0-9]+)$", "\\1", df_ds$sample)
        )
      )
      
      cond <- dplyr::case_when(
        grepl("^SOD1_HEKp_", df_ds$sample) ~ "Parental",
        grepl("^SOD1_KO1\\.2_", df_ds$sample) ~ "Partial KO (1.2)",
        grepl("^SOD1_KO11F6_", df_ds$sample) ~ "Full KO (11F6)",
        TRUE ~ df_ds$sample
      )
      
      cond <- factor(
        cond,
        levels = c(
          "Parental",
          "Partial KO (1.2)",
          "Full KO (11F6)"
        )
      )
      
      dplyr::mutate(
        df_ds,
        condition = cond,
        rep = repnum
      )
    }


    tpms_long_parsed <- reactive({
      
      all <- tpms_long_all()
      
      if (!nrow(all)) return(all)
      
      lapply(
        split(all, all$dataset),
        parse_conditions
      ) |>
        dplyr::bind_rows()
    })

    dat_points <- reactive({
      req(input$gp_datasets, input$gp_gene)
      gene_input <- trimws(input$gp_gene)
      df <- tpms_long_parsed() |>
        dplyr::filter(dataset %in% input$gp_datasets, gene_name == gene_input)
      if (nrow(df) == 0) {
        shiny::showNotification(sprintf("Warning - Gene '%s' was not found in the selected datasets.", gene_input),
                                type = "error", duration = 5)
      }
      df
    })

    dat_summary <- reactive({
      dat_points() |>
        dplyr::group_by(dataset, condition) |>
        dplyr::summarise(n = dplyr::n(), mean = mean(TPM, na.rm = TRUE),
                         sd = stats::sd(TPM, na.rm = TRUE), .groups = "drop")
    })

    build_plot <- function(points, summary, logy = FALSE, title_txt = "") {
      points$dataset  <- factor(points$dataset,  levels = input$gp_datasets)
      summary$dataset <- factor(summary$dataset, levels = input$gp_datasets)
      lab_ds <- ggplot2::labeller(dataset = ggplot2::as_labeller(
        function(v) vapply(v, pretty_label, "", USE.NAMES = FALSE)))
      p <- ggplot2::ggplot() +
        ggplot2::geom_crossbar(data = summary,
                               ggplot2::aes(x = condition, y = mean, ymin = mean, ymax = mean),
                               width = 0.6, linewidth = 1) +
        ggplot2::geom_errorbar(data = summary,
                               ggplot2::aes(x = condition, ymin = mean - sd, ymax = mean + sd),
                               width = 0.3, linewidth = 0.8) +
        ggplot2::geom_point(data = points,
                            ggplot2::aes(x = condition, y = TPM),
                            position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
                            size = 3.5, color = "#005249", na.rm = TRUE, show.legend = FALSE) +
        ggplot2::labs(title = title_txt,  x = NULL, y = "Transcripts Per Million") +
        theme_Marnie() +
        ggplot2::facet_wrap(dplyr::vars(dataset), scales = "free_x", nrow = 1, labeller = lab_ds)
      if (isTRUE(logy)) p <- p + ggplot2::scale_y_continuous(trans = "log10")
      p
    }

    output$gp_plot <- renderPlot({
      pts <- dat_points(); sms <- dat_summary()
      validate(need(nrow(pts) > 0, "No TPM values for this selection."))
      ds_lab <- if (length(input$gp_datasets) == 1) pretty_label(input$gp_datasets)
      else paste0(length(input$gp_datasets), " datasets")
      build_plot(pts, sms, input$gp_logy, paste(input$gp_gene, "-", ds_lab))
    })

    output$gp_table <- DT::renderDataTable({
      DT::datatable(
        dat_points() |>
          dplyr::arrange(dataset, condition, rep) |>
          dplyr::select(
            gene_id,
            gene_name,
            TPM,
            dataset,
            condition,
            rep,
            sample
          ),
        rownames = FALSE,
        options = list(pageLength = 5, scrollX = TRUE)
      )
    }, server = TRUE)

    safe_stem <- function(ds_vec) {
      if (!length(ds_vec)) return("none")
      if (length(ds_vec) <= 3) paste(ds_vec, collapse = "+") else
        paste0(ds_vec[1], "+", length(ds_vec) - 1, "more")
    }

    output$dl_points <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s_replicates.csv",
                                    safe_stem(input$gp_datasets), input$gp_gene),
      content  = function(file) readr::write_csv(dat_points(), file)
    )
    output$dl_summary <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s_summary.csv",
                                    safe_stem(input$gp_datasets), input$gp_gene),
      content  = function(file) readr::write_csv(dat_summary(), file)
    )
    output$dl_plot <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s.svg",
                                    safe_stem(input$gp_datasets), input$gp_gene),
      content = function(file) {
        pts <- isolate(dat_points()); sms <- isolate(dat_summary())
        validate(need(nrow(pts) > 0, "No data available for this gene."))
        ds_lab <- isolate(if (length(input$gp_datasets) == 1)
          pretty_label(input$gp_datasets)
          else paste0(length(input$gp_datasets), " datasets"))
        grDevices::svg(file, width = 11, height = 7, onefile = TRUE)
        print(build_plot(pts, sms, isolate(input$gp_logy),
                         paste(isolate(input$gp_gene), "-", ds_lab)))
        grDevices::dev.off()
      }
    )
    output$dl_plot_png <- shiny::downloadHandler(
      filename = function() sprintf(
        "TPM_%s_%s.png",
        safe_stem(input$gp_datasets), input$gp_gene
      ),
      content = function(file) {

        pts <- isolate(dat_points())
        sms <- isolate(dat_summary())

        validate(need(nrow(pts) > 0, "No data available for this gene."))

        ds_lab <- isolate(
          if (length(input$gp_datasets) == 1) {
            pretty_label(input$gp_datasets)
          } else {
            paste0(length(input$gp_datasets), " datasets")
          }
        )

        gp <- build_plot(
          pts, sms, isolate(input$gp_logy),
          paste(isolate(input$gp_gene), "-", ds_lab)
        )

        ggplot2::ggsave(
          filename = file,
          plot = gp,
          width = 11,
          height = 7,
          dpi = 300
        )
      }
    )

  })
}
