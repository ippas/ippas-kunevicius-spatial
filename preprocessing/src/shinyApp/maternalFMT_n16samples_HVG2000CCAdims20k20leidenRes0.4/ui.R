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

        .control-card,
        .statistics-card {
          border: 1px solid #e2e2e2;
          border-radius: 6px;
          padding: 18px;
          background: #ffffff;
        }

        .control-card h4,
        .statistics-card h4 {
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
          margin-bottom: 3px;
        }

        .control-note {
          margin-top: 7px;
          font-size: 12px;
          color: #999999;
          line-height: 1.35;
        }

        .statistics-note {
          margin-top: 7px;
          font-size: 12px;
          color: #999999;
        }

        .spatial-plot-container {
          width: 100%;
          margin: 0;
          padding: 0;
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
            column(
              width = 3,

              div(
                class = "control-card",

                h4("Gene expression"),

                selectizeInput(
                  inputId = "spatial_gene",
                  label = "Gene",
                  choices = NULL,
                  selected = NULL,
                  options = list(
                    placeholder = "Select gene..."
                  )
                ),

                sliderInput(
                  inputId = "spatial_point_size",
                  label = "Spot size",
                  min = 0.4,
                  max = 4.0,
                  value = 1.3,
                  step = 0.1
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
                    width = 7,
                    actionButton(
                      inputId = "update_indicated_clusters",
                      label = "Update indicated clusters",
                      width = "100%",
                      class = "btn-primary"
                    )
                  ),
                  column(
                    width = 5,
                    actionButton(
                      inputId = "clear_indicated_clusters",
                      label = "Clear indication",
                      width = "100%"
                    )
                  )
                ),

                div(
                  class = "control-note",
                  "All spots use the red expression scale by default. Indicated clusters use the green expression scale."
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
                )
              )
            ),

            column(
              width = 9,

              div(
                class = "spatial-plot-container",
                plotOutput(
                  outputId = "spatial_gene_plot",
                  height = "auto"
                )
              ),

              div(
                class = "statistics-card",

                h4("Selected gene statistics"),

                tableOutput(
                  outputId = "spatial_gene_statistics"
                ),

                div(
                  class = "statistics-note",
                  textOutput(
                    outputId = "spatial_statistics_note",
                    inline = TRUE
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
