# shiny_modules/mod_DEGSEA.R
mod_anchored_gsea_ui <- function(id) {
  source("/home/quentin/01_PROJETS/43_Shiny_bis/R_analysis/anchored_GSEA.R")
  ns <- NS(id)

  page_sidebar(
    sidebar = sidebar(
      accordion(
        id = ns("anchored_gsea"),
        open = FALSE,
        multiple = TRUE,

        accordion_panel(
          "Files *",
          fileInput(ns("bulk_file"),   "Bulk path : *",   accept = c(".csv", ".tsv")),
          fileInput(ns("clinic_file"), "Clinic path : *", accept = c(".csv", ".tsv")),
          shinyFiles::shinyDirButton(ns("output_dir"), "Output directory : *", "Upload")
        ),

        accordion_panel(
          "Filtering on clinic",
          selectInput(
            ns("clinic_filt"), "Filter on clin column :",
            choices = character(0),
            multiple = FALSE, selectize = FALSE, size = 7
          ),
          selectInput(
            ns("kept_modality"), "Kept modality :",
            choices = character(0),
            multiple = TRUE, selectize = FALSE, size = 7
          )
        ),

        accordion_panel(
          "Filtering on gene",
          selectizeInput(
            ns("gene_filt"), "Filter on gene :",
            choices = character(0),
            multiple = FALSE,
            options = list(placeholder = "Type a gene symbol…")
          ),
          numericInput(
            ns("quantile"), "Quantile threshold :",
            value = 0, min = 0, max = 1, step = 0.05
          ),
          selectInput(
            ns("low_or_high"), "Keep low or high expression :",
            choices = c("low", "high"),
            multiple = FALSE, selectize = FALSE, size = 2
          )
        ),

        accordion_panel(
          "Parameters *",
          selectizeInput(
            ns("anchor_gene"), "Select anchor gene :",
            choices = character(0),
            multiple = FALSE,
            options = list(placeholder = "Type a gene symbol…")
          ),
          selectInput(
            ns("adjust_on_cell_fraction"), "Adjust gene expression using ESTIMATE scores :",
            choices = c("immune_score", "stromal_score", "estimate_score"),
            multiple = TRUE, selectize = FALSE, size = 3
          ),
          selectInput(
            ns("GSEA_geneset"), "GSEA geneset : *",
            choices = c("all_pathways", "REACTOME_pathways", "GOBP_pathways", "KEGG_pathways"),
            multiple = FALSE, selectize = FALSE, size = 4
          )
        )
      ),
      actionButton(ns("run_anchored_GSEA"), "Run anchored GSEA")
    ),

    navset_card_tab(
      nav_panel("Logs", card(verbatimTextOutput(ns("logs")))),
      nav_panel("Volcano", card(plotOutput(ns("volcano"), height = 450))),
      nav_panel("Differential Expression", card(DT::DTOutput(ns("top_DE")))),
      nav_panel("GSEA table",  card(DT::DTOutput(ns("top_GSEA")))),
      nav_panel("GESECA table", card(DT::DTOutput(ns("top_GESECA")))),
      nav_panel("ssGSEA table", card(DT::DTOutput(ns("top_ssGSEA"))))
    )
  )
}

mod_anchored_gsea_server <- function(id, roots = c(home = "~")) {
  moduleServer(id, function(input, output, session) {

    clinic_df <- reactive({
      read_delim_auto(input$clinic_file)
    })

    bulk_df <- reactive({
      read_delim_auto(input$bulk_file)
    })

    # update gene list when bulk changes
    observeEvent(input$bulk_file, {
      req(bulk_df())
      genes <- rownames(bulk_df())
      updateSelectizeInput(session, "gene_filt", choices = genes, server = TRUE)
      updateSelectizeInput(session, "anchor_gene", choices = genes, server = TRUE)
    })

    # On met les covariables d'ajustement sous forme de string
    adjust_on_cell_fraction_str <- reactive({
      x <- input$adjust_on_cell_fraction
      if (is.null(x) || length(x) == 0) {
          return(NULL)
      }
      paste(x, collapse = "+")
    })

    # update clinic column lists
    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "clinic_filt", choices = names(clinic_df()), selected = "")
    })

    # shinyFiles directory chooser (IMPORTANT: pass session=)
    shinyFiles::shinyDirChoose(
      input,
      id = "output_dir",
      session = session,
      roots = roots,
      filetypes = c("", "txt", "bigWig", "tsv", "csv", "bw")
    )

    output_dir_path <- reactive({
      req(input$output_dir)
      shinyFiles::parseDirPath(roots, input$output_dir)
    })


    # clinic filter -> modalities
    observeEvent(input$clinic_filt, {
      req(clinic_df(), input$clinic_filt)
      v <- clinic_df()[[input$clinic_filt]]
      mods <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      mods <- mods[!is.na(mods) & mods != ""]
      updateSelectInput(session, "kept_modality", choices = mods, selected = NULL)
    })

    # ---- RUN anchored_GSEA (this must NOT be inside observeEvent) ----
    anchored_gsea_res <- eventReactive(input$run_anchored_GSEA, {
      req(
        clinic_df(), bulk_df(), input$anchor_gene,
        input$GSEA_geneset, output_dir_path()
      )

      withProgress(message = "Running anchored GSEA...", value = 0, {
        incProgress(0.2)
        res <- anchored_GSEA_pipe(
          rnafilt_counts     = bulk_df(),
          clinic_annot       = clinic_df(),
          output_dir         = output_dir_path(),
          adjust_on_ESTIMATE = adjust_on_cell_fraction_str(),
          anchor_gene        = input$anchor_gene,
          pathways_to_use    = input$GSEA_geneset,
          filter_on_clin_col = input$clinic_filt,
          kept_modality      = input$kept_modality,
          filter_by_gene     = input$gene_filt,
          quantile           = input$quantile,
          keep_low_or_high   = input$low_or_high
        )
        incProgress(1)
        res
      })
    })

    # outputs
    output$logs <- renderPrint({
      req(anchored_gsea_res())
      names(anchored_gsea_res())
    })

    output$top_GSEA <- DT::renderDT({
      df <- anchored_gsea_res
      req(df)
      df
    })
  })
}