#!/usr/bin/env Rscript

library(shiny)
library(ggplot2)
library(plotly)

reference_celltype_colors <- c(
  "C-TAC1-KCNQ5" = "#F8766D",
  "C-TAC1-LRP1B" = "#E88526",
  "C-LTMR" = "#D39200",
  "Aβ-HTMR" = "#B79F00",
  "C-OSMR-SST" = "#93AA00",
  "Aδ-COOL" = "#5EB300",
  "C-TAC1-TRPA1" = "#00BA38",
  "Aβ-LTMR-PALM" = "#00BF74",
  "Aδ-TRPV1" = "#00C19F",
  "C-COLD" = "#00BFC4",
  "Aβ-PROPRIO" = "#00B9E3",
  "Aβ-LTMR-SLIT2" = "#00ADFA",
  "C-OSMR-GFRA1_2" = "#619CFF",
  "Aδ-LTMR" = "#DB72FB",
  "Aβ-LTMR-ALDH1A1" = "#F564E3",
  "Aδ-CACNA1E" = "#FF61C3",
  "patchseq" = "#c00000"
)

expr_mat_cpm <- readRDS("data/expr_mat_cpm.rds")
cell_meta <- readRDS("data/cell_meta.rds")
agg_ephys <- readRDS("data/agg_ephys.rds")
agg_gene_cpm <- readRDS("data/agg_gene_cpm_median.rds")
gene_list <- readRDS("data/gene_list.rds")
efeat_list <- readRDS("data/efeat_list.rds")

if (!"labels" %in% colnames(cell_meta)) {
  stop("cell_meta.rds must include a labels column.")
}

ui <- fluidPage(
  titlePanel("Patch-seq gene expression ↔ Patch-clamp electrophysiology feature Explorer"),
  tags$style(HTML("
    .plot-note {
      background: #f5f5f5;
      border: 1px solid #e0e0e0;
      border-radius: 6px;
      padding: 8px 10px;
      margin: 6px 0 10px 0;
      color: #333333;
      font-size: 0.95em;
    }
  ")),
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        "gene",
        "Gene",
        choices = sort(gene_list),
        selected = if ("SCN10A" %in% gene_list) "SCN10A" else sort(gene_list)[1]
      ),
      selectizeInput(
        "efeat",
        "Ephys feature",
        choices = sort(efeat_list),
        selected = if ("Rheobase" %in% efeat_list) "Rheobase" else sort(efeat_list)[1]
      ),
      radioButtons(
        "mode",
        "Scatter plot",
        choices = c("Cell-level" = "cell", "Aggregated by cell-type" = "agg"),
        selected = "agg"
      )
    ),
    mainPanel(
      plotlyOutput("scatter", height = "600px"),
      div(
        class = "plot-note",
        "Selected ephys feature vs. log2(CPM+1) expression.",
        "Colors indicate cell types; trendline shown in black."
      ),
      tags$hr(),
      plotlyOutput("violin", height = "600px"),
      div(
        class = "plot-note",
        "Distribution of log2(CPM+1) expression across cell types. ",
        "Points show individual cells with jitter."
      )
    )
  )
)

server <- function(input, output, session) {
  selected_data <- reactive({
    req(input$gene, input$efeat)

    if (input$mode == "cell") {
      x <- as.numeric(expr_mat_cpm[input$gene, ])
      y <- cell_meta[[input$efeat]]
      data.frame(
        gene = x,
        ephys = y,
        labels = cell_meta$labels,
        check.names = FALSE
      )
    } else {
      if (is.null(agg_gene_cpm)) {
        stop("agg_gene_cpm_median.rds is missing; re-run data/prep_app_data.R with compute_gene_medians = TRUE.")
      }
      x <- as.numeric(agg_gene_cpm[input$gene, ])
      y <- agg_ephys[[input$efeat]]
      data.frame(
        gene = x,
        ephys = y,
        labels = agg_ephys$labels,
        check.names = FALSE
      )
    }
  })

  output$scatter <- renderPlotly({
    df <- selected_data()
    df$gene_log2p1 <- log2(df$gene + 1)
    cor_val <- suppressWarnings(cor(df$gene_log2p1, df$ephys, method = "spearman", use = "pairwise.complete.obs"))
    cor_label <- paste0("Spearman r = ", round(cor_val, 4))
    df$hover <- paste0(
      "label: ", df$labels, "<br>",
      "gene (CPM): ", signif(df$gene, 4), "<br>",
      "efeat: ", signif(df$ephys, 4)
    )

    p <- ggplot(df, aes(x = gene_log2p1, y = ephys, color = labels, text = hover)) +
      geom_point(alpha = 0.8, size = 2) +
      geom_smooth(
        aes(x = gene_log2p1, y = ephys, group = 1),
        inherit.aes = FALSE,
        method = "lm",
        se = FALSE,
        color = "black"
      ) +
      scale_color_manual(values = reference_celltype_colors, na.value = "grey70", name = "Cell types") +
      labs(x = paste0("Gene expr. (log2 CPM+1): ", input$gene), y = input$efeat) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")

    ggplotly(p, tooltip = "text") %>%
      layout(
        annotations = list(
          list(
            x = 1,
            y = 1,
            xref = "paper",
            yref = "paper",
            text = cor_label,
            showarrow = FALSE,
            xanchor = "right",
            yanchor = "top"
          )
        )
      )
  })

  output$violin <- renderPlotly({
    req(input$gene)
    df <- data.frame(
      gene_log2p1 = log2(as.numeric(expr_mat_cpm[input$gene, ]) + 1),
      labels = cell_meta$labels,
      check.names = FALSE
    )

    p <- ggplot(df, aes(x = labels, y = gene_log2p1, color = labels)) +
      geom_violin(fill = NA, linewidth = 0.6, alpha = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
      scale_color_manual(values = reference_celltype_colors, na.value = "grey70", name = "Cell types") +
      labs(x = "Cell types", y = paste0("Gene expr. (log2 CPM+1): ", input$gene)) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    ggplotly(p, tooltip = c("y", "x"))
  })
}

shinyApp(ui, server)
