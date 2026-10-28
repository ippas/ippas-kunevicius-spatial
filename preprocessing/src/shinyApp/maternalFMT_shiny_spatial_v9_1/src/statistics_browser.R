# ==============================================================================
# Differential-expression statistics browser
# ============================================================================== 

statistics_browser_test_catalog <- function() {
  list(
    overall_group = list(
      label = "Donor group main effect",
      aliases = c(
        "Overall_Group_ASD_vs_Neurotypical"
      )
    ),
    overall_sex = list(
      label = "Sex main effect",
      aliases = c(
        "Overall_Sex_Female_vs_Male"
      )
    ),
    interaction = list(
      label = "Group × Sex interaction",
      aliases = c(
        "Interaction"
      )
    ),
    group_male = list(
      label = "ASD vs Neurotypical — Male",
      aliases = c(
        "ASD_vs_Neurotypical_in_Male",
        "ASD_Male_vs_Neurotypical_Male"
      )
    ),
    group_female = list(
      label = "ASD vs Neurotypical — Female",
      aliases = c(
        "ASD_vs_Neurotypical_in_Female",
        "ASD_Female_vs_Neurotypical_Female"
      )
    ),
    sex_neurotypical = list(
      label = "Female vs Male — Neurotypical",
      aliases = c(
        "Neurotypical_Female_vs_Neurotypical_Male",
        "Female_vs_Male_in_Neurotypical"
      )
    ),
    sex_asd = list(
      label = "Female vs Male — ASD",
      aliases = c(
        "ASD_Female_vs_ASD_Male",
        "Female_vs_Male_in_ASD"
      )
    )
  )
}

resolve_statistics_browser_tests <- function(
  edger_results,
  test_definitions
) {
  catalog <- statistics_browser_test_catalog()
  available_ids <- names(edger_results)

  resolved <- lapply(
    names(catalog),
    function(key) {
      entry <- catalog[[key]]
      matching_ids <- entry$aliases[
        entry$aliases %in% available_ids
      ]

      if (length(matching_ids) == 0L) {
        stop(
          "Could not identify edgeR result table for statistics tab '",
          entry$label,
          "'. Expected one of: ",
          paste(entry$aliases, collapse = ", "),
          ". Available tables: ",
          paste(available_ids, collapse = ", "),
          call. = FALSE
        )
      }

      test_id <- matching_ids[[1]]

      comparison <- test_id
      if (
        is.data.frame(test_definitions) &&
          all(c("test_id", "comparison") %in% colnames(test_definitions))
      ) {
        definition_rows <- test_definitions[
          as.character(test_definitions$test_id) == test_id,
          ,
          drop = FALSE
        ]
        if (nrow(definition_rows) == 1L) {
          comparison_value <- as.character(
            definition_rows$comparison[[1]]
          )
          if (!is.na(comparison_value) && nzchar(comparison_value)) {
            comparison <- comparison_value
          }
        }
      }

      list(
        key = key,
        label = entry$label,
        test_id = test_id,
        comparison = comparison
      )
    }
  )

  names(resolved) <- names(catalog)
  resolved
}

pretty_statistics_group_label <- function(x) {
  x <- gsub("_", " ", as.character(x), fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

format_mean_sd_statistics <- function(
  mean_value,
  sd_value,
  digits = 1L
) {
  mean_value <- as.numeric(mean_value)
  sd_value <- as.numeric(sd_value)

  output <- rep(NA_character_, length(mean_value))

  both_finite <- is.finite(mean_value) & is.finite(sd_value)
  mean_only <- is.finite(mean_value) & !is.finite(sd_value)

  output[both_finite] <- paste0(
    formatC(mean_value[both_finite], format = "f", digits = digits),
    " ± ",
    formatC(sd_value[both_finite], format = "f", digits = digits)
  )

  output[mean_only] <- formatC(
    mean_value[mean_only],
    format = "f",
    digits = digits
  )

  output
}

add_statistics_summary_columns <- function(
  table,
  show_percent_positive = FALSE,
  show_tmm_cpm = FALSE
) {
  output <- table

  if (isTRUE(show_percent_positive)) {
    mean_columns <- grep(
      "^mean_percent_positive_spots_",
      colnames(table),
      value = TRUE
    )

    for (mean_column in mean_columns) {
      suffix <- sub(
        "^mean_percent_positive_spots_",
        "",
        mean_column
      )
      sd_column <- paste0(
        "sd_percent_positive_spots_",
        suffix
      )

      if (sd_column %in% colnames(table)) {
        display_name <- paste0(
          "Positive spots (%) — ",
          pretty_statistics_group_label(suffix)
        )
        output[[display_name]] <- format_mean_sd_statistics(
          mean_value = table[[mean_column]],
          sd_value = table[[sd_column]],
          digits = 1L
        )
      }
    }
  }

  if (isTRUE(show_tmm_cpm)) {
    mean_columns <- grep(
      "^mean_TMM_CPM_",
      colnames(table),
      value = TRUE
    )

    for (mean_column in mean_columns) {
      suffix <- sub(
        "^mean_TMM_CPM_",
        "",
        mean_column
      )
      sd_column <- paste0(
        "sd_TMM_CPM_",
        suffix
      )

      if (sd_column %in% colnames(table)) {
        display_name <- paste0(
          "TMM-CPM — ",
          pretty_statistics_group_label(suffix)
        )
        output[[display_name]] <- format_mean_sd_statistics(
          mean_value = table[[mean_column]],
          sd_value = table[[sd_column]],
          digits = 2L
        )
      }
    }
  }

  output
}

prepare_statistics_browser_table <- function(
  result_table,
  selected_clusters,
  abs_log2fc_threshold,
  significance_metric,
  significance_threshold,
  selected_columns,
  show_percent_positive,
  show_tmm_cpm
) {
  table <- as.data.frame(
    result_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "cluster_id",
    "logFC",
    "PValue",
    "FDR"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Statistics table is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  selected_clusters <- as.character(selected_clusters)

  table <- table[
    as.character(table$cluster_id) %in% selected_clusters,
    ,
    drop = FALSE
  ]

  if (is.finite(abs_log2fc_threshold) && abs_log2fc_threshold > 0) {
    keep_logfc <- is.finite(table$logFC) &
      abs(table$logFC) >= abs_log2fc_threshold
    table <- table[keep_logfc, , drop = FALSE]
  }

  significance_metric <- match.arg(
    significance_metric,
    c("FDR", "PValue")
  )

  significance_values <- as.numeric(
    table[[significance_metric]]
  )
  keep_significance <- is.finite(significance_values) &
    significance_values <= significance_threshold
  table <- table[keep_significance, , drop = FALSE]

  table <- add_statistics_summary_columns(
    table = table,
    show_percent_positive = show_percent_positive,
    show_tmm_cpm = show_tmm_cpm
  )

  base_column_map <- c(
    cluster_id = "Cluster",
    ensembl_gene_id = "Ensembl gene ID",
    gene = "Gene",
    chromosome = "Chromosome",
    logFC = "log2FC",
    logCPM = "logCPM",
    F = "F",
    PValue = "PValue",
    FDR = "FDR",
    regulation = "Regulation"
  )

  selected_base_columns <- selected_columns[
    selected_columns %in% colnames(table)
  ]

  summary_columns <- colnames(table)[
    grepl(
      "^(Positive spots \\(%\\)|TMM-CPM) — ",
      colnames(table)
    )
  ]

  display_columns <- c(
    selected_base_columns,
    summary_columns
  )

  display_columns <- unique(
    display_columns[
      display_columns %in% colnames(table)
    ]
  )

  if (length(display_columns) == 0L) {
    stop(
      "No display columns were selected.",
      call. = FALSE
    )
  }

  display <- table[
    ,
    display_columns,
    drop = FALSE
  ]

  if ("cluster_id" %in% colnames(display)) {
    display$cluster_id <- paste0(
      "C",
      as.character(display$cluster_id)
    )
  }

  rename_columns <- intersect(
    names(base_column_map),
    colnames(display)
  )

  for (column_name in rename_columns) {
    colnames(display)[
      colnames(display) == column_name
    ] <- base_column_map[[column_name]]
  }

  list(
    display = display,
    n_filtered = nrow(table),
    n_total = nrow(result_table)
  )
}

statistics_browser_ui <- function(id) {
  ns <- shiny::NS(id)

  test_catalog <- statistics_browser_test_catalog()

  table_tabs <- lapply(
    names(test_catalog),
    function(key) {
      shiny::tabPanel(
        title = test_catalog[[key]]$label,
        shiny::div(
          class = "statistics-comparison-note",
          shiny::textOutput(
            outputId = ns(paste0("comparison_note_", key)),
            inline = TRUE
          )
        ),
        shiny::div(
          class = "statistics-result-count",
          shiny::textOutput(
            outputId = ns(paste0("result_count_", key)),
            inline = TRUE
          )
        ),
        DT::DTOutput(
          outputId = ns(paste0("table_", key))
        )
      )
    }
  )

  shiny::fluidRow(
    shiny::column(
      width = 3,
      shiny::div(
        class = "statistics-browser-sidebar",

        shiny::h4("Filters"),

        shiny::numericInput(
          inputId = ns("abs_log2fc"),
          label = "Minimum |log2FC|",
          value = 0.5,
          min = 0,
          step = 0.1
        ),
        shiny::div(
          class = "parameter-note",
          "Minimum absolute effect size."
        ),

        shiny::fluidRow(
          shiny::column(
            width = 6,
            shiny::selectInput(
              inputId = ns("significance_metric"),
              label = "Significance",
              choices = c(
                "FDR" = "FDR",
                "PValue" = "PValue"
              ),
              selected = "FDR"
            )
          ),
          shiny::column(
            width = 6,
            shiny::numericInput(
              inputId = ns("significance_threshold"),
              label = "Threshold",
              value = 0.05,
              min = 0,
              max = 1,
              step = 0.01
            )
          )
        ),
        shiny::div(
          class = "parameter-note",
          "Maximum FDR or nominal P-value."
        ),

        shiny::hr(),

        shiny::h4("Clusters"),
        shiny::div(
          class = "statistics-cluster-checkboxes",
          shiny::checkboxGroupInput(
            inputId = ns("clusters"),
            label = NULL,
            choices = character(0),
            selected = character(0)
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 6,
            shiny::actionButton(
              inputId = ns("select_all_clusters"),
              label = "Select all",
              width = "100%"
            )
          ),
          shiny::column(
            width = 6,
            shiny::actionButton(
              inputId = ns("clear_clusters"),
              label = "Clear",
              width = "100%"
            )
          )
        ),
        shiny::div(
          class = "parameter-note",
          "Only selected clusters are shown."
        ),

        shiny::hr(),

        shiny::h4("Displayed columns"),
        shiny::div(
          class = "statistics-column-checkboxes",
          shiny::checkboxGroupInput(
            inputId = ns("columns"),
            label = NULL,
            choices = c(
              "Cluster ID" = "cluster_id",
              "Ensembl gene ID" = "ensembl_gene_id",
              "Gene" = "gene",
              "Chromosome" = "chromosome",
              "log2FC" = "logFC",
              "logCPM" = "logCPM",
              "F" = "F",
              "PValue" = "PValue",
              "FDR" = "FDR",
              "Regulation" = "regulation"
            ),
            selected = c(
              "cluster_id",
              "gene",
              "F",
              "PValue",
              "FDR",
              "regulation"
            )
          )
        ),

        shiny::checkboxInput(
          inputId = ns("show_percent_positive"),
          label = "Positive spots (%) — mean ± SD",
          value = FALSE
        ),
        shiny::checkboxInput(
          inputId = ns("show_tmm_cpm"),
          label = "TMM-CPM — mean ± SD",
          value = FALSE
        ),
        shiny::div(
          class = "parameter-note",
          "Group summaries are added only when selected."
        )
      )
    ),

    shiny::column(
      width = 9,
      shiny::div(
        class = "statistics-browser-main",
        do.call(
          shiny::tabsetPanel,
          c(
            list(
              id = ns("test_tab"),
              type = "tabs"
            ),
            table_tabs
          )
        ),
        shiny::div(
          class = "statistics-browser-footnote",
          shiny::tags$b("TMM-CPM: "),
          "counts per million after edgeR TMM normalization for library composition. ",
          shiny::tags$b("Positive spots: "),
          "percentage of spots with detected expression; group values are shown as mean ± SD."
        )
      )
    )
  )
}

statistics_browser_server <- function(
  id,
  edger_results,
  test_definitions,
  cluster_ids,
  cluster_names
) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      resolved_tests <- resolve_statistics_browser_tests(
        edger_results = edger_results,
        test_definitions = test_definitions
      )

      cluster_ids <- sort_cluster_ids(cluster_ids)

      cluster_labels <- vapply(
        cluster_ids,
        function(cluster_id) {
          anatomical_name <- cluster_names[[as.character(cluster_id)]]
          if (
            is.null(anatomical_name) ||
              is.na(anatomical_name) ||
              !nzchar(anatomical_name)
          ) {
            paste0("C", cluster_id)
          } else {
            paste0(
              "C",
              cluster_id,
              " — ",
              anatomical_name
            )
          }
        },
        character(1)
      )

      cluster_choices <- stats::setNames(
        cluster_ids,
        cluster_labels
      )

      shiny::updateCheckboxGroupInput(
        session = session,
        inputId = "clusters",
        choices = cluster_choices,
        selected = cluster_ids
      )

      shiny::observeEvent(
        input$select_all_clusters,
        {
          shiny::updateCheckboxGroupInput(
            session = session,
            inputId = "clusters",
            choices = cluster_choices,
            selected = cluster_ids
          )
        }
      )

      shiny::observeEvent(
        input$clear_clusters,
        {
          shiny::updateCheckboxGroupInput(
            session = session,
            inputId = "clusters",
            choices = cluster_choices,
            selected = character(0)
          )
        }
      )

      for (key in names(resolved_tests)) {
        local({
          current_key <- key
          current_test <- resolved_tests[[current_key]]

          output[[paste0("comparison_note_", current_key)]] <-
            shiny::renderText({
              paste0(
                current_test$comparison,
                " · test_id: ",
                current_test$test_id
              )
            })

          prepared_table <- shiny::reactive({
            shiny::req(input$significance_metric)
            shiny::req(input$significance_threshold)
            shiny::req(input$abs_log2fc)

            shiny::validate(
              shiny::need(
                length(input$clusters) > 0L,
                "Select at least one cluster."
              ),
              shiny::need(
                length(input$columns) > 0L ||
                  isTRUE(input$show_percent_positive) ||
                  isTRUE(input$show_tmm_cpm),
                "Select at least one column to display."
              ),
              shiny::need(
                is.finite(input$abs_log2fc) &&
                  input$abs_log2fc >= 0,
                "Minimum |log2FC| must be 0 or greater."
              ),
              shiny::need(
                is.finite(input$significance_threshold) &&
                  input$significance_threshold >= 0 &&
                  input$significance_threshold <= 1,
                "Significance threshold must be between 0 and 1."
              )
            )

            prepare_statistics_browser_table(
              result_table = edger_results[[current_test$test_id]],
              selected_clusters = input$clusters,
              abs_log2fc_threshold = input$abs_log2fc,
              significance_metric = input$significance_metric,
              significance_threshold = input$significance_threshold,
              selected_columns = input$columns,
              show_percent_positive = input$show_percent_positive,
              show_tmm_cpm = input$show_tmm_cpm
            )
          })

          output[[paste0("result_count_", current_key)]] <-
            shiny::renderText({
              result <- prepared_table()
              paste0(
                "Filtered results: ",
                format(result$n_filtered, big.mark = ","),
                " / ",
                format(result$n_total, big.mark = ",")
              )
            })

          output[[paste0("table_", current_key)]] <- DT::renderDT(
            {
              result <- prepared_table()
              display <- result$display

              shiny::validate(
                shiny::need(
                  nrow(display) > 0L,
                  "No results meet the current filters."
                )
              )

              dt <- DT::datatable(
                display,
                rownames = FALSE,
                filter = "top",
                class = "stripe hover compact",
                options = list(
                  pageLength = 25,
                  lengthMenu = c(10, 25, 50, 100),
                  lengthChange = TRUE,
                  scrollX = TRUE,
                  autoWidth = TRUE,
                  deferRender = TRUE,
                  order = list(),
                  dom = "lftip"
                )
              )

              if ("log2FC" %in% colnames(display)) {
                dt <- DT::formatRound(
                  dt,
                  columns = "log2FC",
                  digits = 3
                )
              }

              if ("logCPM" %in% colnames(display)) {
                dt <- DT::formatRound(
                  dt,
                  columns = "logCPM",
                  digits = 3
                )
              }

              if ("F" %in% colnames(display)) {
                dt <- DT::formatRound(
                  dt,
                  columns = "F",
                  digits = 3
                )
              }

              probability_columns <- intersect(
                c("PValue", "FDR"),
                colnames(display)
              )

              if (length(probability_columns) > 0L) {
                dt <- DT::formatSignif(
                  dt,
                  columns = probability_columns,
                  digits = 3
                )
              }

              dt
            },
            server = TRUE
          )
        })
      }

      invisible(resolved_tests)
    }
  )
}
