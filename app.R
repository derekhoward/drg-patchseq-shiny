#!/usr/bin/env Rscript

library(shiny)
library(ggplot2)
library(plotly)
library(Seurat)

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
integrated_obj <- readRDS("data/drg_integrated_slim.RDS")

if (!"umap" %in% names(integrated_obj@reductions)) {
  stop("drg_integrated.RDS must include a UMAP reduction named 'umap'.")
}

integrated_meta <- integrated_obj@meta.data
if (!all(c("labels", "labels.p", "dataset") %in% colnames(integrated_meta))) {
  missing_cols <- setdiff(c("labels", "labels.p", "dataset"), colnames(integrated_meta))
  stop(paste0("drg_integrated.RDS is missing required metadata columns: ", paste(missing_cols, collapse = ", "), "."))
}

umap_embeddings <- integrated_obj@reductions$umap@cell.embeddings
umap_df <- data.frame(
  cell_id = rownames(umap_embeddings),
  UMAP_1 = umap_embeddings[, 1],
  UMAP_2 = umap_embeddings[, 2],
  labels = integrated_meta$labels,
  "labels.p" = integrated_meta[["labels.p"]],
  dataset = integrated_meta$dataset,
  check.names = FALSE
)
umap_df$dataset_display <- ifelse(umap_df$dataset == "patch", "patchseq", "snRNAseq")
umap_df$is_patch <- umap_df$dataset == "patch"
umap_df$hover <- paste0(
  "dataset: ", umap_df$dataset_display, "<br>",
  "cell-type: ", umap_df$labels
)
umap_limits <- list(
  x = range(umap_df$UMAP_1, na.rm = TRUE),
  y = range(umap_df$UMAP_2, na.rm = TRUE)
)

if (!"labels" %in% colnames(cell_meta)) {
  stop("cell_meta.rds must include a labels column.")
}

ui <- navbarPage(
  title = "Patch-seq Explorer",
  header = tags$style(HTML("
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
  tabPanel(
    "Home",
    fluidPage(
      h2("Patch-seq gene expression ↔ electrophysiology feature explorer"),
      p("Use the Gene Search tab to explore gene expression relationships with electrophysiology features and integrated UMAPs."),
      tags$hr(),
      h3("The Data"),
      p("Unlike standard transcriptomics, this study utilizes Patch-seq to link gene expression directly to cellular function. By combining electrophysiological characterization with single-cell sequencing, we have mapped the transcriptomic profile of functionally identified CMis."),
      h3("Our Objective"),
      p("To characterize the molecular identity of human dermal nociceptors and provide a mechanistic basis for understanding neuropathic pain. This platform integrates our core Patch-seq findings with supporting snRNA-seq data to facilitate the discovery of novel therapeutic targets.")
    )
  ),
  tabPanel(
    "Gene Search",
    titlePanel("Patch-seq gene expression ↔ electrophysiology feature explorer"),
    sidebarLayout(
      sidebarPanel(
        selectizeInput(
          "gene",
          "Gene",
          choices = NULL
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
          selected = "cell"
        ),
        tags$hr(),
        plotlyOutput("umap", height = "350px"),
        div(
          class = "plot-note",
          "Integrated pig DRG UMAP showing both snRNAseq and Patch-seq datasets."
        ),
        tags$hr(),
        plotlyOutput("feature_umap", height = "420px")
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
)

server <- function(input, output, session) {
  observeEvent(TRUE, {
    updateSelectizeInput(
      session,
      "gene",
      choices = sort(gene_list),
      selected = if ("SCN10A" %in% gene_list) "SCN10A" else sort(gene_list)[1],
      server = TRUE
    )
  }, once = TRUE)

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

    plot_obj <- ggplotly(p, tooltip = "text") %>%
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
        ),
        hoverlabel = list(align = "left")
      )
    plot_obj <- plotly::style(plot_obj, hoverinfo = "text")
    plotly::toWebGL(plot_obj)
  })

  output$violin <- renderPlotly({
    req(input$gene)
    df <- data.frame(
      gene_log2p1 = log2(as.numeric(expr_mat_cpm[input$gene, ]) + 1),
      labels = cell_meta$labels,
      check.names = FALSE
    )

    p <- ggplot(df, aes(x = labels, y = gene_log2p1, color = labels)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.6, alpha = 0.4) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
      scale_color_manual(values = reference_celltype_colors, na.value = "grey70", name = "Cell types") +
      labs(x = "Cell types", y = paste0("Gene expr. (log2 CPM+1): ", input$gene)) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    plotly::style(ggplotly(p, tooltip = "none"), hoverinfo = "skip")
  })

  output$umap <- renderPlotly({
    snrna_df <- umap_df[!umap_df$is_patch, , drop = FALSE]
    patch_df <- umap_df[umap_df$is_patch, , drop = FALSE]

    p <- ggplot() +
      geom_point(
        data = snrna_df,
        aes(x = UMAP_1, y = UMAP_2, color = `labels.p`, text = hover),
        size = 0.6,
        alpha = 0.65
      ) +
      geom_point(
        data = patch_df,
        aes(x = UMAP_1, y = UMAP_2, fill = `labels.p`, text = hover),
        shape = 21,
        size = 1.5,
        stroke = 0.4,
        color = "black",
        alpha = 0.9
      ) +
      scale_color_manual(values = reference_celltype_colors, na.value = "grey70", name = "Predicted cell types") +
      scale_fill_manual(values = reference_celltype_colors, na.value = "grey70", guide = "none") +
      labs(x = "UMAP 1", y = "UMAP 2") +
      coord_equal(xlim = umap_limits$x, ylim = umap_limits$y) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9)
      )

    plotly::toWebGL(
      ggplotly(p, tooltip = "text") %>%
      layout(
        margin = list(l = 35, r = 10, t = 10, b = 30),
        hoverlabel = list(align = "left")
      )
    )
  })

  output$feature_umap <- renderPlotly({
    req(input$gene)

    assay_name <- "integrated"
    validate(
      need(
        assay_name %in% names(integrated_obj@assays),
        "Assay 'integrated' not found in drg_integrated.RDS."
      )
    )
    validate(
      need(
        input$gene %in% rownames(integrated_obj@assays[[assay_name]]),
        paste0("Gene ", input$gene, " not found in assay ", assay_name, ".")
      )
    )

    prev_assay <- Seurat::DefaultAssay(integrated_obj)
    Seurat::DefaultAssay(integrated_obj) <- assay_name
    on.exit(Seurat::DefaultAssay(integrated_obj) <- prev_assay, add = TRUE)

    feature_plot <- Seurat::FeaturePlot(
      integrated_obj,
      features = input$gene,
      reduction = "umap",
      cols = c("grey90", "#d73027"),
      pt.size = 0.6,
      order = TRUE,
      combine = FALSE
    )
    p <- feature_plot[[1]]

    plot_data <- p$data
    plot_data$cell_type <- integrated_meta[rownames(plot_data), "labels"]
    plot_data$hover_text <- paste0("cell-type: ", plot_data$cell_type)
    p$data <- plot_data

    p <- p +
      aes(text = hover_text) +
      labs(title = input$gene, x = "UMAP_1", y = "UMAP_2") +
      coord_equal(xlim = umap_limits$x, ylim = umap_limits$y) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 6)),
        plot.title.position = "plot"
      )

    plotly::toWebGL(
      ggplotly(p, tooltip = "text") %>%
      layout(
        margin = list(l = 35, r = 10, t = 40, b = 30),
        hoverlabel = list(align = "left")
      )
    )
  })
}

shinyApp(ui, server)
