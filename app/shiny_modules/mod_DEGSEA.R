# shiny_modules/mod_DEGSEA.R
mod_degsea_ui <- function(id) {
  ns <- NS(id)

  page_sidebar(
    sidebar = sidebar(
      accordion(
        id = ns("degsea_acc"),
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
          selectInput(
            ns("DESeq_covar"), "DESeq covariates :",
            choices = character(0),
            multiple = TRUE, selectize = FALSE, size = 7
          ),
          selectInput(
            ns("DESeq_contrast"), "DESeq contrast : *",
            choices = character(0),
            multiple = FALSE, selectize = FALSE, size = 7
          ),
          selectInput(
            ns("DESeq_control"), "DESeq control : *",
            choices = character(0),
            multiple = TRUE, selectize = FALSE, size = 5
          ),
          selectInput(
            ns("DESeq_test"), "DESeq test : *",
            choices = character(0),
            multiple = TRUE, selectize = FALSE, size = 5
          ),
          selectInput(
            ns("GSEA_geneset"), "GSEA geneset : *",
            choices = c("all_pathways", "REACTOME_pathways", "GOBP_pathways", 
                        "KEGG_pathways", "QoL_pathways", "transcript_factor_target",
                        "boyault_sets"),
            multiple = FALSE, selectize = FALSE, size = 4
          )
        )
      ),
      actionButton(ns("run_DEGSEA"), "Run DEGSEA")
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

mod_degsea_server <- function(id, roots = c(home = "~")) {
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
    })

    # update clinic column lists
    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "DESeq_contrast", choices = names(clinic_df()), selected = "")
      updateSelectInput(session, "DESeq_covar",    choices = names(clinic_df()), selected = NULL)
      updateSelectInput(session, "clinic_filt",    choices = names(clinic_df()), selected = "")
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

    # contrast -> control/test modalities
    observeEvent(input$DESeq_contrast, {
      req(clinic_df(), input$DESeq_contrast)
      v <- clinic_df()[[input$DESeq_contrast]]
      mods <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      mods <- mods[!is.na(mods) & mods != ""]
      updateSelectInput(session, "DESeq_control", choices = mods, selected = "")
      updateSelectInput(session, "DESeq_test",    choices = mods, selected = "")
    })

    # clinic filter -> modalities
    observeEvent(input$clinic_filt, {
      req(clinic_df(), input$clinic_filt)
      v <- clinic_df()[[input$clinic_filt]]
      mods <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      mods <- mods[!is.na(mods) & mods != ""]
      updateSelectInput(session, "kept_modality", choices = mods, selected = NULL)
    })

    # ---- RUN DEGSEA (this must NOT be inside observeEvent) ----
    degsea_res <- eventReactive(input$run_DEGSEA, {
      req(
        clinic_df(), bulk_df(),
        input$DESeq_contrast, input$DESeq_control, input$DESeq_test,
        input$GSEA_geneset, output_dir_path()
      )

      withProgress(message = "Running DEGSEA...", value = 0, {
        incProgress(0.2)
        res <- DEGSEA_pipe(
          rnafilt_counts     = bulk_df(),
          clinic_annot       = clinic_df(),
          contrast           = input$DESeq_contrast,
          control            = input$DESeq_control,
          test               = input$DESeq_test,
          pathways_to_use    = input$GSEA_geneset,
          output_dir         = output_dir_path(),
          covariates         = input$DESeq_covar,
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
      req(degsea_res())
      names(degsea_res())
    })

    output$volcano <- renderPlot({
      p <- degsea_res()$plots$volcano
      req(p)
      print(p)
    })

    output$top_DE <- DT::renderDT({
      df <- degsea_res()$tables$top_DE
      req(df)
      df
    })

    output$top_GSEA <- DT::renderDT({
      df <- degsea_res()$tables$top_GSEA
      req(df)
      df
    })

    output$top_GESECA <- DT::renderDT({
      df <- degsea_res()$tables$top_GESECA
      req(df)
      df
    })

    output$top_ssGSEA <- DT::renderDT({
      df <- degsea_res()$tables$top_ssGSEA
      req(df)
      df
    })
  })
}