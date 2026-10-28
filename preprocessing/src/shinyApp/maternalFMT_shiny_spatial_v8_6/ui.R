# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        body {
          background-color: #ffffff;
        }

        .app-container {
          width: calc(100% - 24px);
          max-width: none;
          margin: 0 auto;
          padding: 18px 12px 24px 12px;
        }

        .app-header {
          margin-bottom: 25px;
        }

        .app-title {
          font-size: 32px;
          font-weight: 600;
          margin: 0;
        }

        .app-subtitle {
          margin-top: 6px;
          font-size: 15px;
          color: #5f6368;
        }

        .app-description {
          margin-top: 4px;
          font-size: 12px;
          color: #a0a0a0;
        }

        .nav-tabs {
          font-size: 16px;
        }

        .nav-tabs > li > a {
          padding: 12px 20px;
        }

        .tab-content {
          padding-top: 22px;
        }

        .spatial-subtabs {
          margin-top: 0;
        }

        .spatial-subtabs > .nav-tabs {
          font-size: 14px;
          margin-bottom: 16px;
        }

        .spatial-subtabs > .nav-tabs > li > a {
          padding: 8px 15px;
        }

        .control-card,
        .statistics-card,
        .barplots-card {
          border: 1px solid #e2e2e2;
          border-radius: 6px;
          padding: 18px;
          background: #ffffff;
        }

        .control-card h4,
        .statistics-card h4,
        .barplots-card h4 {
          margin-top: 0;
          font-weight: 600;
        }

        .cluster-checkboxes .shiny-options-group {
          display: grid;
          grid-template-columns: 1fr 1fr;
          column-gap: 12px;
        }

        .cluster-checkboxes .checkbox {
          margin-top: 3px;
          margin-bottom: 5px;
        }

        .cluster-checkboxes .checkbox label {
          font-size: 12px;
          line-height: 1.15;
          font-weight: 400;
        }

        .control-note,
        .statistics-note {
          margin-top: 7px;
          font-size: 12px;
          color: #999999;
          line-height: 1.4;
        }

        .statistics-explanation {
          margin-top: 12px;
          font-size: 12px;
          color: #777777;
          line-height: 1.45;
        }

        .spatial-expression-note {
          padding-top: 8px;
          font-size: 12px;
          color: #777777;
          line-height: 1.4;
        }

        .spatial-plot-container {
          width: 100%;
          margin: 0;
          padding: 0;
        }

        .compact-numeric .form-group {
          margin-bottom: 0;
        }

        .tab-placeholder {
          padding: 40px 10px;
          color: #777777;
        }
      ")
    )
  ),

  div(
    class = "app-container",

    div(
      class = "app-header",

      div(
        class = "app-title",
        "maternalFMT ASD — Spatial Transcriptomics"
      ),

      div(
        class = "app-subtitle",
        "16 samples · Clustering criteria: 2,000 highly variable genes (HVG) · ",
        "CCA integration · 20 dimensions · k = 20 nearest neighbors · ",
        "Leiden resolution = 0.4"
      ),

      div(
        class = "app-description",
        "HVG — genes used for clustering · ",
        "CCA — sample integration · ",
        "dimensions — components used for neighbor graph · ",
        "k — number of nearest neighbors · ",
        "resolution — clustering granularity"
      )
    ),

    tabsetPanel(
      id = "main_tab",

      tabPanel(
        title = "Spatial visualization",
        icon = icon("brain"),

        div(
          class = "tab-content",

          fluidRow(
            # ------------------------------------------------------------------
            # Persistent controls: unchanged when switching right-side subtabs
            # ------------------------------------------------------------------
            column(
              width = 3,

              div(
                class = "control-card",

                h4("Gene expression"),

                fluidRow(
                  column(
                    width = 8,
                    selectizeInput(
                      inputId = "spatial_gene",
                      label = "Gene",
                      choices = NULL,
                      selected = NULL,
                      options = list(
                        placeholder = "Select gene..."
                      )
                    )
                  ),
                  column(
                    width = 4,
                    div(
                      class = "compact-numeric",
                      numericInput(
                        inputId = "spatial_point_size",
                        label = "Spot size",
                        value = 1.3,
                        min = 0.1,
                        max = 6,
                        step = 0.1
                      )
                    )
                  )
                ),

                div(
                  class = "control-note",
                  "Spot size controls dot size on the plot."
                ),

                fluidRow(
                  column(
                    width = 6,
                    numericInput(
                      inputId = "spatial_upper_quantile",
                      label = "Upper colour quantile",
                      value = 0.99,
                      min = 0.5,
                      max = 1,
                      step = 0.01
                    )
                  )
                ),

                div(
                  class = "control-note",
                  "Colour cap: 0.99 clips the top 1% of expression values so extreme spots do not dominate the scale."
                ),

                hr(),

                h4("Indicate clusters"),

                div(
                  class = "cluster-checkboxes",
                  checkboxGroupInput(
                    inputId = "spatial_clusters",
                    label = NULL,
                    choices = character(0),
                    selected = character(0)
                  )
                ),

                fluidRow(
                  column(
                    width = 4,
                    actionButton(
                      inputId = "select_all_clusters",
                      label = "Select all",
                      width = "100%"
                    )
                  ),
                  column(
                    width = 4,
                    actionButton(
                      inputId = "update_indicated_clusters",
                      label = "Update",
                      width = "100%",
                      class = "btn-primary"
                    )
                  ),
                  column(
                    width = 4,
                    actionButton(
                      inputId = "clear_indicated_clusters",
                      label = "Clear",
                      width = "100%"
                    )
                  )
                ),

                div(
                  class = "control-note",
                  "Indicated clusters use the red expression scale; all remaining clusters use the green expression scale."
                ),

                fluidRow(
                  column(
                    width = 6,
                    textAreaInput(
                      inputId = "red_palette_input",
                      label = "Red scale (indicated)",
                      value = "#D9D9D9, #FEE5D9, #FCAE91, #FB6A4A, #DE2D26, #A50F15, #67000D",
                      rows = 3,
                      width = "100%"
                    )
                  ),
                  column(
                    width = 6,
                    textAreaInput(
                      inputId = "green_palette_input",
                      label = "Green scale (remaining)",
                      value = "#D9D9D9, #E5F5E0, #A1D99B, #74C476, #31A354, #006D2C, #00441B",
                      rows = 3,
                      width = "100%"
                    )
                  )
                ),

                div(
                  class = "control-note",
                  "Comma-separated hex colours; first colour = low expression, last colour = high expression."
                ),

                hr(),

                h4("Significance"),

                selectInput(
                  inputId = "spatial_effect",
                  label = "Effect",
                  choices = c(
                    "Donor group main effect" = "group",
                    "Sex main effect" = "sex",
                    "Group × Sex interaction" = "interaction"
                  ),
                  selected = "group"
                ),

                fluidRow(
                  column(
                    width = 6,
                    numericInput(
                      inputId = "spatial_fdr",
                      label = "FDR",
                      value = 0.05,
                      min = 0,
                      max = 1,
                      step = 0.01
                    )
                  ),
                  column(
                    width = 6,
                    numericInput(
                      inputId = "spatial_abs_log2fc",
                      label = "|log2FC|",
                      value = 0.5,
                      min = 0,
                      step = 0.1
                    )
                  )
                ),

                actionButton(
                  inputId = "indicate_significant_clusters",
                  label = "Indicate significant clusters",
                  width = "100%"
                ),

                div(
                  class = "control-note",
                  "FDR and |log2FC| define which clusters are automatically selected as significant."
                )
              )
            ),

            # ------------------------------------------------------------------
            # Right display area: only this area changes between subtabs
            # ------------------------------------------------------------------
            column(
              width = 9,

              div(
                class = "spatial-subtabs",

                tabsetPanel(
                  id = "spatial_subtab",
                  type = "tabs",

                  tabPanel(
                    title = "Spatial plot",
                    icon = icon("image"),

                    div(
                      style = "text-align:right; margin-bottom:6px;",
                      downloadButton(
                        outputId = "download_spatial_pdf",
                        label = "Download PDF",
                        class = "btn-default"
                      )
                    ),

                    div(
                      class = "spatial-plot-container",
                      plotOutput(
                        outputId = "spatial_gene_plot",
                        height = "auto"
                      )
                    )
                  ),

                  tabPanel(
                    title = "Cluster expression",
                    icon = icon("layer-group"),

                    div(
                      style = "text-align:right; margin-bottom:6px;",
                      downloadButton(
                        outputId = "download_cluster_spatial_pdf",
                        label = "Download PDF",
                        class = "btn-default"
                      )
                    ),

                    div(
                      class = "spatial-plot-container",
                      plotOutput(
                        outputId = "cluster_spatial_gene_plot",
                        height = "auto"
                      )
                    )
                  ),

                  tabPanel(
                    title = "Statistics info",
                    icon = icon("table"),

                    div(
                      class = "statistics-card",

                      uiOutput(
                        outputId = "spatial_statistics_title"
                      ),

                      div(
                        class = "statistics-note",
                        textOutput(
                          outputId = "spatial_statistics_note",
                          inline = TRUE
                        )
                      ),

                      br(),

                      DT::DTOutput(
                        outputId = "spatial_gene_statistics"
                      ),

                      div(
                        class = "statistics-explanation",
                        tags$b("Table values: "),
                        "positive spots are shown as mean ± SD percentage of spots with detected expression within each comparison group. ",
                        "TMM-CPM is shown as mean ± SD counts per million after edgeR TMM normalization for between-sample library composition."
                      )
                    )
                  ),

                  tabPanel(
                    title = "Barplots",
                    icon = icon("chart-column"),

                    div(
                      class = "barplots-card",
                      h4("Barplots"),
                      p("Gene-expression summary barplots will be displayed here.")
                    )
                  )
                )
              )
            )
          )
        )
      ),

      tabPanel(
        title = "Heatmaps",
        icon = icon("th"),

        div(
          class = "tab-placeholder",
          h3("Heatmaps"),
          p("Visualization of expression patterns for selected genes.")
        )
      ),

      tabPanel(
        title = "Statistics",
        icon = icon("chart-bar"),

        div(
          class = "tab-placeholder",
          h3("Statistics"),
          p("Differential gene expression results.")
        )
      )
    )
  )
)
