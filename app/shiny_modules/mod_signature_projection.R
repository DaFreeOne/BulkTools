# shiny_modules/mod_signature_projection.R
library(jsonlite)

mod_signature_proj_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML(sprintf("
      #%s {
        width: 100%%;
        height: 100%%;
        display: flex;
        align-items: center;
        justify-content: center;
        overflow: hidden;
      }

      #%s img {
        max-width: 100%%;
        max-height: 100%%;
        width: auto;
        height: auto;
        object-fit: fill;    
        display: block;
      }
    ", ns("projection"), ns("projection")))),

    page_sidebar(
      sidebar = sidebar(
        accordion(
          id = ns("sign_proj_acc"),
          open = FALSE,
          multiple = TRUE,

          accordion_panel(
            "Files *",
            fileInput(ns("bulk_file"),   "Reference bulk : *",   accept = c(".csv", ".tsv")),
            fileInput(ns("clinic_file"), "Reference clinic : *", accept = c(".csv", ".tsv")),
            fileInput(ns("proj_sample"), "Sample to project : *", accept = c(".csv", ".tsv")),
            shinyFiles::shinyDirButton(ns("output_dir"), "Output directory : *", "Upload")
          ),
          accordion_panel(
            "Filtering on clinic",
            selectInput(
              ns("clinic_filt_proj"), "Filter on clin column :",
              choices = character(0),
              multiple = FALSE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("kept_modality_proj"), "Kept modality :",
              choices = character(0),
              multiple = TRUE, selectize = FALSE, size = 7
            )
          ),
          accordion_panel(
            "Filtering on gene",
            selectizeInput(
              ns("gene_filt_proj"), "Filter on gene :",
              choices = character(0),
              multiple = FALSE,
              options = list(placeholder = "Type a gene symbol…")
            ),
            numericInput(
              ns("quantile_proj"), "Quantile threshold :",
              value = 0, min = 0, max = 1, step = 0.05
            ),
            selectInput(
              ns("low_or_high_proj"), "Keep low or high expression :",
              choices = c("low", "high"),
              multiple = FALSE, selectize = FALSE, size = 2
            )
          ),
          accordion_panel(
            "Select signature *",
            selectizeInput(
              ns("therapy"), "Therapy :",
              choices = character(0),
              multiple = FALSE,
              options = list(placeholder = "Type a treatment (no spaces)…")
            ),
            selectizeInput(
              ns("signatures"), "Signature :",
              choices = character(0),
              multiple = FALSE,
              options = list(placeholder = "Type a signature name (no spaces)…")
            )
          ),
          accordion_panel(
            "Parameters *",
            selectInput(
              ns("contrast"), "Contrast :",
              choices = character(0),
              multiple = TRUE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("responders"), "Responders :",
              choices = character(0),
              multiple = TRUE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("non_responders"), "Non-Responders :",
              choices = character(0),
              multiple = TRUE, selectize = FALSE, size = 7
            )
          ),
          accordion_panel(
            "KM plot parameters ",
            selectInput(
              ns("plot_km_or_not"), "Plot KM plot ? :",
              choices = c("Ja", "Nein"),
              multiple = FALSE, selectize = FALSE, size = 2
            ),
            selectInput(
              ns("event_realization_col"), "Event realization column :",
              choices = character(0),
              multiple = FALSE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("survival_time_col"), "Survival time column :",
              choices = character(0),
              multiple = FALSE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("group_quantile"), "Quantile stratification :",
              choices = c("median", "tertile", "quartile"),
              multiple = FALSE, selectize = FALSE, size = 3
            ),
          )
        ),
        actionButton(ns("run_signature_proj"), "Run signature projection")
      ),
      navset_card_tab(
        nav_panel("Logs", card(verbatimTextOutput(ns("logs")))),
        nav_panel("Projection", card(imageOutput(ns("projection"), height = "575px"))),    
        nav_panel("Survival analysis", card(imageOutput(ns("KM_plot"), height = 500))),
        nav_panel("Signature assessment", 
                  layout_columns(
                    col_width = c(6, 6, 6),
                    card(imageOutput(ns("signature_assessment_boxplot"), width = "100%", height = "300px")),
                    card(imageOutput(ns("signature_assessment_ROC"), width = "100%", height = "300px"))
                    # card(imageOutput(ns("signature_assessment_conf_mat"), width = "100%", height = "300px"))
                  )
        )
      )    
    )
  )  
}


mod_signature_proj_server <- function(id, roots = c(home = "~")) {
  moduleServer(id, function(input, output, session) {
    library(reticulate)
    Sys.unsetenv("RETICULATE_PYTHON")
    use_condaenv("myShiny", required = TRUE)
    source_python("py/py_plots.py")


    clinic_df <- reactive({
      read_delim_auto(input$clinic_file)
    })

    bulk_df <- reactive({
      read_delim_auto(input$bulk_file)
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

    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "contrast", choices = names(clinic_df()), selected = "")
      updateSelectInput(session, "survival_time_col", choices = names(clinic_df()),  selected = "") 
      updateSelectInput(session, "event_realization_col", choices = names(clinic_df()),  selected = "") 
      updateSelectInput(session, "clinic_filt_proj",    choices = names(clinic_df()), selected = "")
    })

    # update gene list when bulk changes
    observeEvent(input$bulk_file, {
      req(bulk_df())
      genes <- rownames(bulk_df())
      updateSelectizeInput(session, "gene_filt_proj", choices = genes, server = TRUE)
    })

    # clinic filter -> modalities
    observeEvent(input$clinic_filt_proj, {
      req(clinic_df(), input$clinic_filt_proj)
      v <- clinic_df()[[input$clinic_filt_proj]]
      mods <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      mods <- mods[!is.na(mods) & mods != ""]
      updateSelectInput(session, "kept_modality_proj", choices = mods, selected = NULL)
    })

    observeEvent(input$contrast, {
      req(clinic_df(), input$contrast)
      v <- clinic_df()[[input$contrast]]
      mods <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      mods <- mods[!is.na(mods) & mods != ""]
      updateSelectInput(session, "responders", choices = mods, selected = "")
      updateSelectInput(session, "non_responders",    choices = mods, selected = "")
    })

    app_dir <- normalizePath(getwd())
    signatures_path = file.path(app_dir, "REF_DATA", "signatures", "signatures.json")
    signatures_json = fromJSON(signatures_path)
    updateSelectizeInput(session, "therapy", choices = names(signatures_json), server = TRUE)

    observeEvent(input$therapy, {
      req(input$therapy, signatures_json)
      signatures_list <- names(signatures_json[[input$therapy]])
      updateSelectizeInput(session, "signatures", choices = signatures_list, server = TRUE)
    })

    sel_signature <- reactive({
      req(input$therapy, input$signatures, signatures_json)
      signatures_json[[input$therapy]][[input$signatures]]
    })

    projection_res <- eventReactive(input$run_signature_proj, {
      req(
        clinic_df(), bulk_df(), sel_signature(),
        input$proj_sample$datapath, output_dir_path()
      )

      withProgress(message = "Projecting signature score...", value = 0, {
        incProgress(0.2)
        output_dir = file.path(path.expand(output_dir_path()), input$therapy, input$signatures)
        res <- signature_proj_pipe(
          rnafilt_counts         = bulk_df(),
          clinic_annot           = clinic_df(),
          output_dir             = output_dir,
          therapy_used           = input$therapy,
          signature_name         = input$signatures,
          signature_to_use       = sel_signature()$geneset,
          sample_to_project_path = input$proj_sample$datapath[1],
          contrast               = input$contrast,
          resp_var               = input$responders,
          non_resp_var           = input$non_responders, 
          survival_time_col      = input$survival_time_col,
          event_realization_col  = input$event_realization_col,
          group_quantile         = input$group_quantile,
          filter_on_clin_col     = input$clinic_filt_proj,
          kept_modality          = input$kept_modality_proj,
          filter_by_gene         = input$gene_filt_proj,
          keep_low_or_high       = input$low_or_high_proj,
          quantile_thr           = input$quantile_proj
        )
        incProgress(1)
        res
      })
    })

    output$logs <- renderPrint({
      req(projection_res())
      names(projection_res())
    })

    output$projection <- renderImage({
      width  <- session$clientData$output_myplot_width
      height <- session$clientData$output_myplot_height
      p <- projection_res()$proj_plot$save_path
      req(p, file.exists(p))
      list(
        src = p,
        contentType = "image/png",
        alt = "Signature projection plot"
      )
    }, deleteFile = FALSE)

    output$KM_plot <- renderImage({
      width  <- session$clientData$output_myplot_width
      height <- session$clientData$output_myplot_height
      p <- projection_res()$KM_plot$save_path
      req(p, file.exists(p))
      list(
        src = p,
        width = width,
        height = height,
        contentType = "image/png",
        alt = "Signature survival KM"
      )
    }, deleteFile = FALSE)

    output$signature_assessment_ROC <- renderImage({
      width  <- session$clientData$output_myplot_width
      height <- session$clientData$output_myplot_height
      roc_p <- projection_res()$roc_plot$roc_save_path
      req(roc_p, file.exists(roc_p))
      list(
        src = roc_p,
        width = width,
        height = height,
        contentType = "image/png",
        alt = "Signature evaluation ROC"
      )

    }, deleteFile = FALSE)
    output$signature_assessment_boxplot <- renderImage({
      width  <- session$clientData$output_myplot_width
      height <- session$clientData$output_myplot_height
      box_p <- projection_res()$box_plot$save_path
      req(box_p, file.exists(box_p))
      list(
        src = box_p,
        width = width,
        height = height,
        contentType = "image/png",
        alt = "Signature evaluation boxplot"
      )
    }, deleteFile = FALSE)

  })
}
