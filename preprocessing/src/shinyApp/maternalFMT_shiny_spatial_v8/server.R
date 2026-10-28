# ==============================================================================
# Server
# ==============================================================================

server <- function(input, output, session) {
  seurat_object <- app_data$seurat_object
  edger_results <- app_data$edger_results
  edger_test_definitions <- app_data$edger_test_definitions

  sample_layout <- get_sample_layout_metadata(
    seurat_object = seurat_object,
    sample_id_column = sample_id_column,
    group_column = group_column,
    sex_column = sex_column
  )

  available_samples <- as.character(
    sample_layout$sample_ID
  )

  available_clusters <- get_available_cluster_ids(
    seurat_object = seurat_object,
    cluster_column = cluster_column
  )

  cluster_choices <- stats::setNames(
    available_clusters,
    paste0("C", available_clusters)
  )

  default_gene <- if ("Drd2" %in% app_data$genes) {
    "Drd2"
  } else if ("Glul" %in% app_data$genes) {
    "Glul"
  } else {
    app_data$genes[[1]]
  }

  updateSelectizeInput(
    session = session,
    inputId = "spatial_gene",
    choices = app_data$genes,
    selected = default_gene,
    server = TRUE
  )

  # All clusters are indicated at application start.
  updateCheckboxGroupInput(
    session = session,
    inputId = "spatial_clusters",
    choices = cluster_choices,
    selected = available_clusters
  )

  effect_test_ids <- resolve_shiny_effect_test_ids(
    test_definitions = edger_test_definitions,
    edger_results = edger_results
  )

  message("Selected edgeR tests for Shiny:")
  message("  donor group: ", effect_test_ids[["group"]])
  message("  sex: ", effect_test_ids[["sex"]])
  message("  interaction: ", effect_test_ids[["interaction"]])

  section_contexts <- prepare_all_spatial_sections(
    seurat_object = seurat_object,
    sample_ids = available_samples,
    sample_id_column = sample_id_column,
    cluster_column = cluster_column,
    image_scale = image_scale
  )

  indicated_clusters <- reactiveVal(
    available_clusters
  )

  observeEvent(
    input$select_all_clusters,
    {
      updateCheckboxGroupInput(
        session = session,
        inputId = "spatial_clusters",
        choices = cluster_choices,
        selected = available_clusters
      )
      indicated_clusters(available_clusters)
    }
  )

  observeEvent(
    input$update_indicated_clusters,
    {
      indicated_clusters(
        as.character(input$spatial_clusters)
      )
    }
  )

  observeEvent(
    input$clear_indicated_clusters,
    {
      updateCheckboxGroupInput(
        session = session,
        inputId = "spatial_clusters",
        choices = cluster_choices,
        selected = character(0)
      )
      indicated_clusters(character(0))
    }
  )

  observeEvent(
    input$indicate_significant_clusters,
    {
      req(input$spatial_gene)
      req(input$spatial_effect)

      significant_clusters <- get_significant_clusters_for_gene(
        edger_results = edger_results,
        effect_test_ids = effect_test_ids,
        effect = input$spatial_effect,
        gene_symbol = input$spatial_gene,
        fdr_threshold = input$spatial_fdr,
        abs_log2fc_threshold = input$spatial_abs_log2fc
      )

      message(
        "Significant clusters | gene = ",
        input$spatial_gene,
        " | effect = ",
        input$spatial_effect,
        " | FDR < ",
        input$spatial_fdr,
        " | |log2FC| >= ",
        input$spatial_abs_log2fc,
        " | clusters = ",
        if (length(significant_clusters) == 0L) {
          "<none>"
        } else {
          paste(
            paste0("C", significant_clusters),
            collapse = ", "
          )
        }
      )

      updateCheckboxGroupInput(
        session = session,
        inputId = "spatial_clusters",
        choices = cluster_choices,
        selected = significant_clusters
      )

      indicated_clusters(significant_clusters)
    }
  )

  spatial_plot_data <- reactive({
    req(input$spatial_gene)

    add_gene_expression_to_all_sections(
      section_contexts = section_contexts,
      seurat_object = seurat_object,
      gene = input$spatial_gene,
      expression_assay = expression_assay,
      expression_layer_prefix = expression_layer_prefix,
      counts_layer_prefix = counts_layer_prefix,
      sample_order = available_samples
    )
  })

  spatial_plot_height <- reactive({
    plot_width <- session$clientData$output_spatial_gene_plot_width
    n_rows <- max(
      as.integer(
        table(sample_layout$display_group)
      )
    )

    if (is.null(plot_width) || !is.finite(plot_width) || plot_width <= 0) {
      return(1750)
    }

    panel_width <- plot_width / 4

    as.integer(
      max(
        1050,
        round(
          panel_width * n_rows +
            42 * max(n_rows - 1L, 0L) +
            235
        )
      )
    )
  })

  output$spatial_gene_plot <- renderPlot(
    {
      req(input$spatial_gene)
      req(input$spatial_point_size)

      validate(
        need(
          is.finite(input$spatial_point_size) &&
            input$spatial_point_size > 0,
          "Spot size must be greater than 0."
        )
      )

      create_gene_spatial_plot_all_sections(
        plot_data = spatial_plot_data(),
        section_contexts = section_contexts,
        sample_layout = sample_layout,
        gene = input$spatial_gene,
        indicated_clusters = indicated_clusters(),
        red_palette = red_expression_palette,
        green_palette = green_expression_palette,
        point_size = input$spatial_point_size
      )
    },
    height = function() {
      spatial_plot_height()
    },
    res = 110
  )

  selected_gene_statistics <- reactive({
    req(input$spatial_gene)
    req(input$spatial_effect)

    get_gene_edger_statistics(
      edger_results = edger_results,
      edger_test_definitions = edger_test_definitions,
      effect_test_ids = effect_test_ids,
      effect = input$spatial_effect,
      gene_symbol = input$spatial_gene,
      cluster_names = cluster_names
    )
  })

  output$spatial_statistics_title <- renderUI({
    req(input$spatial_gene)

    h4(
      paste0(
        "Selected gene statistics — ",
        input$spatial_gene
      )
    )
  })

  output$spatial_gene_statistics <- DT::renderDT(
    {
      statistics <- selected_gene_statistics()

      validate(
        need(
          nrow(statistics$display) > 0L,
          "No edgeR results were found for the selected gene."
        )
      )

      dt <- DT::datatable(
        statistics$display,
        rownames = FALSE,
        filter = "top",
        class = "stripe hover compact",
        options = list(
          pageLength = 16,
          lengthChange = FALSE,
          scrollX = TRUE,
          autoWidth = TRUE,
          order = list(),
          dom = "tip"
        )
      )

      dt <- DT::formatRound(
        dt,
        columns = "log2FC",
        digits = 3
      )

      dt <- DT::formatSignif(
        dt,
        columns = c("PValue", "FDR"),
        digits = 3
      )

      dt
    },
    server = FALSE
  )

  output$spatial_statistics_note <- renderText({
    statistics <- selected_gene_statistics()

    paste0(
      statistics$comparison,
      " · test_id: ",
      statistics$result_id
    )
  })
}
