gseaExploreUI <- function(id) {
  ns <- shiny::NS(id)
  
  tagList(
    
    checkboxGroupInput(
      ns("datasets"),
      "Comparisons",
      choices = character(0)
    ),
    
    selectInput(
      ns("ontology"),
      "Gene set source",
      choices = c("All", "BP", "CC", "MF", "KEGG", "Reactome", "TF_GTRD"),
      selected = "All",
      multiple = TRUE
    ),
    
    checkboxInput(
      ns("sig_only"),
      "Show only FDR < 0.05",
      value = TRUE
    ),
    
    radioButtons(
      ns("direction"),
      "Direction",
      choices = c(
        "Both" = "both",
        "Positive NES" = "up",
        "Negative NES" = "down"
      ),
      selected = "both",
      inline = TRUE
    ),
    
    numericInput(
      ns("min_sig"),
      "Term must be significant in at least",
      value = 1,
      min = 1,
      step = 1
    ),
    
    numericInput(
      ns("max_plot_terms"),
      "Maximum terms to plot",
      value = 15,
      min = 5,
      max = 500,
      step = 5
    ),
    selectInput(
      ns("theme"),
      "Predefined biological category",
      choices = c(
        "None - use custom search" = "custom",
        "Mitochondria" = "mito",
        "Cellular stress" = "stress",
        "Oxidative stress" = "oxidative_stress",
        "ER stress / unfolded protein response" = "ER_stress",
        "Autophagy / lysosome" = "autophagy",
        "Lipid metabolism" = "lipid",
        "Cytoskeleton / axon / transport" = "cytoskeleton",
        "Iron metabolism" = "iron"
      ),
      selected = "custom"
    ),
    
    textInput(
      ns("keyword"),
      "Search GSEA terms",
      placeholder = "e.g. mitochondrial translation, axon, DNA repair"
    ),
    
    checkboxInput(
      ns("regex"),
      "Use advanced pattern matching",
      value = FALSE
    ),
    
    helpText(
      "Optional: use | for OR searches, e.g. mitochondria|respiratory chain. Leave unchecked for a normal phrase search."
    )
  )
}

gseaExploreMainUI <- function(id) {
  ns <- shiny::NS(id)
  
  tagList(
    
    uiOutput(ns("n_terms")),
    
    br(),
    
    shinycssloaders::withSpinner(
      uiOutput(ns("gsea_plot_ui")),
      type = 4,
      color = "#005249"
    ),
    
    br(),
    
    h4("Filtered GSEA Results"),
    
    DT::DTOutput(ns("results")),
    
    div(
      class = "text-center mt-2",
      downloadButton(
        ns("dl_gsea_plot_svg"),
        "Download plot (SVG)"
      )
    ),
    
    div(
      class = "text-center mt-2",
      downloadButton(
        ns("dl_results"),
        "Download filtered results (CSV)"
      )
    ),
    
    br()
  )
}


gseaExploreServer <- function(id, pkg = utils::packageName()) {
  
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("stringr", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("forcats", quietly = TRUE)
    requireNamespace("DT", quietly = TRUE)
    
    `%||%` <- function(x, y) {
      if (is.null(x) || length(x) == 0) y else x
    }
    
    
    # =========================================================
    # 1. Load GSEA long-format table once
    # =========================================================
    
    gsea_tbl <- reactiveVal(NULL)
    
    observeEvent(TRUE, {
      
      # CHANGE THIS LINE to whatever cloud loader you use
      x <- get_gsea_explore_data_cloud_cached()
      
      x <- as.data.frame(x)
      
      # Standardise names if necessary
      if ("padj" %in% names(x) && !"p.adjust" %in% names(x)) {
        x$p.adjust <- x$padj
      }
      
      gsea_tbl(x)
      
    }, once = TRUE)
    
    
    # =========================================================
    # 2. Pretty dataset labels
    # =========================================================
    
    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) character(0)
    )
    
    pretty_label_vec <- function(ids) {
      
      ids <- as.character(ids)
      
      out <- unname(pretty_map[ids])
      
      bad <- is.na(out) | !nzchar(out)
      out[bad] <- ids[bad]
      
      out
    }
    
    
    # =========================================================
    # 3. Populate comparison choices
    # =========================================================
    
    observe({
      
      df <- gsea_tbl()
      req(df)
      
      ids <- sort(unique(df$dataset))
      labs <- pretty_label_vec(ids)
      
      updateCheckboxGroupInput(
        session,
        "datasets",
        choices = stats::setNames(ids, labs),
        selected = intersect(
          input$datasets %||% character(0),
          ids
        )
      )
    })
    
    
    # =========================================================
    # 4. Biological theme dictionary
    # =========================================================
    
    pathway_sets <- list(
      
      mito =
        paste(
          c(
            "mitochondr",
            "mitophagy",
            "respiratory",
            "electron transport",
            "oxidative phosphorylation",
            "OXPHOS",
            "ATP synthase",
            "TCA cycle",
            "citric acid cycle"
          ),
          collapse = "|"
        ),
      
      stress =
        paste(
          c(
            "stress response",
            "cellular response to stress",
            "oxidative stress",
            "ER stress",
            "endoplasmic reticulum stress",
            "unfolded protein",
            "heat shock",
            "DNA damage",
            "hypoxia",
            "reactive oxygen"
          ),
          collapse = "|"
        ),
      
      oxidative_stress =
        paste(
          c(
            "oxidative stress",
            "reactive oxygen",
            "ROS",
            "redox",
            "antioxidant",
            "glutathione",
            "peroxide"
          ),
          collapse = "|"
        ),
      
      ER_stress =
        paste(
          c(
            "endoplasmic reticulum",
            "ER stress",
            "unfolded protein response",
            "UPR",
            "protein folding",
            "ERAD"
          ),
          collapse = "|"
        ),
      
      autophagy =
        paste(
          c(
            "autophagy",
            "autophagic",
            "autophagosome",
            "lysosome",
            "lysosomal",
            "mitophagy",
            "macroautophagy",
            "chaperone mediated autophagy"
          ),
          collapse = "|"
        ),
      
      lipid =
        paste(
          c(
            "lipid",
            "fatty acid",
            "cholesterol",
            "sterol",
            "phospholipid",
            "sphingolipid",
            "ceramide",
            "beta oxidation",
            "fatty acid oxidation",
            "cardiolipin",
            "ferroptosis"
          ),
          collapse = "|"
        ),
      
      cytoskeleton =
        paste(
          c(
            "cytoskeleton",
            "microtubule",
            "actin",
            "growth cone",
            "neurite",
            "axon",
            "axonal",
            "dendrite",
            "cell projection",
            "axonogenesis",
            "intracellular transport",
            "organelle transport",
            "kinesin",
            "dynein"
          ),
          collapse = "|"
        ),
      
      iron =
        paste(
          c(
            "iron",
            "iron-sulfur",
            "iron sulfur",
            "ferritin",
            "ferropt",
            "heme",
            "transferrin"
          ),
          collapse = "|"
        )
    )
    
    
    # =========================================================
    # 5. Resolve theme + user keyword into search pattern
    # =========================================================
    search_pattern <- reactive({
      
      theme <- input$theme %||% "custom"
      keyword <- trimws(input$keyword %||% "")
      
      patterns <- character(0)
      
      # Predefined theme
      if (
        theme != "custom" &&
        theme %in% names(pathway_sets)
      ) {
        patterns <- c(patterns, pathway_sets[[theme]])
      }
      
      # User-entered keyword
      if (nzchar(keyword)) {
        
        if (isTRUE(input$regex)) {
          # Advanced regex search
          user_pattern <- keyword
        } else {
          # Normal literal text search
          user_pattern <- stringr::str_replace_all(
            keyword,
            "([.\\\\+*?\\[\\]^$(){}=!<>|:\\-])",
            "\\\\\\1"
          )
        }
        
        patterns <- c(patterns, user_pattern)
      }
      
      if (!length(patterns)) {
        return(NULL)
      }
      
      paste(patterns, collapse = "|")
    })
    
    
    # =========================================================
    # 6. Main filtering
    # =========================================================
    
    filtered_gsea <- reactive({
      
      df <- gsea_tbl()
      req(df)
      
      # Comparisons
      if (length(input$datasets)) {
        df <- df |>
          dplyr::filter(
            .data$dataset %in% input$datasets
          )
      }
      
      # Ontology / database
      if (
        length(input$ontology) &&
        !"All" %in% input$ontology
      ) {
        df <- df |>
          dplyr::filter(
            .data$ontology %in% input$ontology
          )
      }
      
      # Direction
      direction <- input$direction %||% "both"
      
      if (direction == "up") {
        df <- df |>
          dplyr::filter(.data$NES > 0)
      }
      
      if (direction == "down") {
        df <- df |>
          dplyr::filter(.data$NES < 0)
      }
      
      # Theme / keyword search
      pattern <- search_pattern()
      
      if (!is.null(pattern) && nzchar(pattern)) {
        
        df <- df |>
          dplyr::filter(
            stringr::str_detect(
              .data$Description,
              stringr::regex(
                pattern,
                ignore_case = TRUE
              )
            )
          )
      }
      
      # Avoid duplicated terms within a comparison
      df <- df |>
        dplyr::distinct(
          .data$dataset,
          .data$ontology,
          .data$ID,
          .data$Description,
          .keep_all = TRUE
        )
      
      df
    })
    
    
    # =========================================================
    # 7. Determine which pathways to retain
    #
    # This reproduces your:
    # keep_if_sig_in_at_least = n
    # =========================================================
    pathway_keep <- reactive({
      
      df <- filtered_gsea()
      
      min_sig <- input$min_sig %||% 1
      
      if (!nrow(df)) {
        return(
          tibble::tibble(
            ontology = character(),
            ID = character(),
            Description = character(),
            n_present = integer(),
            n_sig = integer()
          )
        )
      }
      
      df |>
        dplyr::group_by(
          .data$ontology,
          .data$ID,
          .data$Description
        ) |>
        dplyr::summarise(
          n_present = dplyr::n_distinct(.data$dataset),
          
          n_sig = sum(
            !is.na(.data$p.adjust) &
              .data$p.adjust < 0.05,
            na.rm = TRUE
          ),
          
          .groups = "drop"
        ) |>
        dplyr::filter(
          .data$n_sig >= min_sig
        )
    })
    
    # =========================================================
    # 8. Final displayed data
    # =========================================================
    display_gsea <- reactive({
      
      df <- filtered_gsea()
      keep <- pathway_keep()
      
      # Add display columns BEFORE potentially returning an empty table
      df <- df |>
        dplyr::mutate(
          comparison = pretty_label_vec(.data$dataset),
          
          sig = !is.na(.data$p.adjust) &
            .data$p.adjust < 0.05,
          
          sig_label = ifelse(
            .data$sig,
            "FDR < 0.05",
            "Not significant"
          ),
          
          size_value = ifelse(
            is.na(.data$p.adjust),
            NA_real_,
            -log10(.data$p.adjust)
          )
        )
      
      # No matching pathways
      if (!nrow(df) || !nrow(keep)) {
        return(df[0, ])
      }
      
      df <- df |>
        dplyr::semi_join(
          keep,
          by = c(
            "ontology",
            "ID",
            "Description"
          )
        )
      
      # Optional strict FDR filter
      if (isTRUE(input$sig_only)) {
        
        df <- df |>
          dplyr::filter(
            !is.na(.data$p.adjust),
            .data$p.adjust < 0.05
          )
      }
      
      df
    })
    
    
    plot_gsea <- reactive({
      
      df <- display_gsea()
      
      validate(
        need(nrow(df) > 0, "No GSEA terms match the current filters.")
      )
      
      max_terms <- input$max_plot_terms %||% 30
      
      # Rank terms by best FDR across selected datasets
      top_terms <- df %>%
        dplyr::group_by(
          ontology,
          ID,
          Description
        ) %>%
        dplyr::summarise(
          best_padj = min(p.adjust, na.rm = TRUE),
          max_abs_NES = max(abs(NES), na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(best_padj, dplyr::desc(max_abs_NES)) %>%
        dplyr::slice_head(n = max_terms)
      
      df %>%
        dplyr::semi_join(
          top_terms,
          by = c("ontology", "ID", "Description")
        )
    })
    # =========================================================
    # 9. Number of matching terms
    # =========================================================
    
    output$n_terms <- renderUI({
      
      all_df <- display_gsea()
      plot_df <- plot_gsea()
      
      n_all <- all_df |>
        dplyr::distinct(
          .data$ontology,
          .data$ID,
          .data$Description
        ) |>
        nrow()
      
      n_plot <- plot_df |>
        dplyr::distinct(
          .data$ontology,
          .data$ID,
          .data$Description
        ) |>
        nrow()
      
      tags$strong(
        sprintf(
          "%d matching GSEA terms; plotting top %d",
          n_all,
          n_plot
        )
      )
    })
    
    # =========================================================
    # 10. Plot
    # =========================================================
    comparison_cols <- c(
      "11F6 full KO vs parental" = "#FE6100",
      "1.2 partial KO vs parental" = "#648FFF",
      "11F6 full KO vs 1.2 partial KO" = "#FFB000"
    )
    
    gsea_plot <- reactive({
      
      df <- plot_gsea()
      
      validate(
        need(nrow(df) > 0, "No GSEA terms match the current filters.")
      )
      
      df <- df |>
        dplyr::mutate(
          Description_wrapped =
            stringr::str_wrap(
              .data$Description,
              width = 48
            )
        )
      
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = .data$NES,
          y = forcats::fct_reorder(
            .data$Description_wrapped,
            .data$NES,
            .fun = mean
          ),
          colour = .data$comparison,
          size = .data$size_value,
          shape = .data$sig_label
        )
      ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dashed"
        ) +
        ggplot2::geom_point(
          position = ggplot2::position_dodge(width = 0.65),
          stroke = 1
        ) +
        ggplot2::scale_shape_manual(
          values = c(
            "FDR < 0.05" = 16,
            "Not significant" = 1
          )
        ) +
        
        ggplot2::scale_colour_manual(
          values = comparison_cols,
          drop = FALSE
        ) +
        ggplot2::labs(
          x = "Normalised Enrichment Score",
          y = NULL,
          colour = NULL,
          shape = NULL,
          size = "-log10 adjusted p-value"
        ) +
        ggplot2::theme_bw(base_size = 26) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = 20)
        )
    })
    
    
    output$gsea_plot <- renderPlot({
      gsea_plot()
    })
    
    
    # =========================================================
    # 11. Results table
    # =========================================================
    
    output$results <- DT::renderDT({
      
      df <- display_gsea()
      
      out <- df |>
        dplyr::transmute(
          Comparison = .data$comparison,
          Source = .data$ontology,
          ID = .data$ID,
          Term = .data$Description,
          NES = round(.data$NES, 3),
          FDR = signif(.data$p.adjust, 3)
        )
      
      DT::datatable(
        out,
        rownames = FALSE,
        filter = "top",
        selection = "multiple",
        options = list(
          pageLength = 25,
          scrollX = TRUE
        )
      )
    })
    
    
    output$gsea_plot_ui <- renderUI({
      
      df <- plot_gsea()
      
      n_terms <- df |>
        dplyr::distinct(
          .data$ontology,
          .data$ID,
          .data$Description
        ) |>
        nrow()
      
      # Approx. 40 px per term plus space for axes/legend
      height_px <- max(
        500,
        180 + n_terms * 40
      )
      
      plotOutput(
        ns("gsea_plot"),
        height = paste0(height_px, "px")
      )
    })
    
    # =========================================================
    # 12. Download filtered data
    # =========================================================
    
    output$dl_results <- downloadHandler(
      
      filename = function() {
        paste0(
          "GSEA_explorer_",
          Sys.Date(),
          ".csv"
        )
      },
      
      content = function(file) {
        
        utils::write.csv(
          display_gsea(),
          file,
          row.names = FALSE
        )
      }
    )
    
    output$dl_gsea_plot_svg <- downloadHandler(
      
      filename = function() {
        paste0(
          "GSEA_filtered_plot_",
          Sys.Date(),
          ".svg"
        )
      },
      
      content = function(file) {
        
        df <- plot_gsea()
        
        n_terms <- df |>
          dplyr::distinct(
            .data$ontology,
            .data$ID,
            .data$Description
          ) |>
          nrow()
        
        # Match exported height to number of plotted terms
        height_in <- max(
          5,
          2 + n_terms * 0.35
        )
        
        svglite::svglite(
          file,
          width = 12,
          height = height_in
        )
        
        print(gsea_plot())
        
        grDevices::dev.off()
      }
    )
    
  })
}

