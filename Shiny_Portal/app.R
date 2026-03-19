#!/usr/bin/env Rscript

current_app_dir <- function() {
  ofiles <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else frame$ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles) > 0) {
    return(dirname(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

app_dir <- current_app_dir()
repo_dir <- normalizePath(file.path(app_dir, ".."), winslash = "/", mustWork = FALSE)

options(cosmx.disable_cli = TRUE)
source(file.path(repo_dir, "CosMx_pipeline.R"), local = TRUE)

suppressPackageStartupMessages({
  library(bslib)
  library(DT)
  library(shiny)
})

pretty_step_label <- function(step_id) {
  labels <- c(
    "01_load" = "01 Load data",
    "02_qc" = "02 Quality control",
    "03_norm" = "03 Normalization",
    "04_dimred" = "04 Dimensionality reduction",
    "05_cluster" = "05 Clustering",
    "06_markers" = "06 Marker analysis",
    "07_annotate" = "07 Annotation",
    "08_visualize" = "08 Visualisation",
    "09_spatial" = "09 Spatial network analysis",
    "10_cci" = "10 CCI analysis",
    "11_batch" = "11 Harmony batch correction",
    "12_bcell" = "12 B-cell microenvironment",
    "merge" = "Merge sample objects"
  )
  labels[[step_id]] %||% step_id
}

tail_text <- function(path, n = 200) {
  if (is.null(path) || !file.exists(path)) {
    return("No log file is available yet.")
  }
  paste(tail(readLines(path, warn = FALSE), n), collapse = "\n")
}

read_table_preview <- function(path, n_max = 1000) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE, n_max = n_max))
  }
  if (ext %in% c("tsv", "txt")) {
    return(readr::read_tsv(path, show_col_types = FALSE, n_max = n_max))
  }
  tibble(note = "Preview is only available for CSV/TSV/TXT files.")
}

list_latest_run_files <- function(output_dir) {
  run_root <- file.path(output_dir, "pipeline_runs")
  if (!dir.exists(run_root)) {
    return(list(latest_run = NULL, sample_results = NULL, merged_results = NULL))
  }
  
  run_dirs <- list.dirs(run_root, recursive = FALSE, full.names = TRUE)
  if (length(run_dirs) == 0) {
    return(list(latest_run = NULL, sample_results = NULL, merged_results = NULL))
  }
  
  latest_run <- run_dirs[which.max(file.info(run_dirs)$mtime)]
  list(
    latest_run = latest_run,
    sample_results = file.path(latest_run, "sample_pipeline_results.csv"),
    merged_results = file.path(latest_run, "merged_pipeline_results.csv")
  )
}

list_target_directories <- function(output_dir) {
  targets <- c()
  labels <- c()
  
  sample_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[basename(sample_dirs) != "pipeline_runs" & basename(sample_dirs) != "Merged" & basename(sample_dirs) != "Paired"]
  if (length(sample_dirs) > 0) {
    targets <- c(targets, sample_dirs)
    labels <- c(labels, basename(sample_dirs))
  }
  
  merged_root <- file.path(output_dir, "Merged")
  if (dir.exists(merged_root)) {
    merged_dirs <- list.dirs(merged_root, recursive = FALSE, full.names = TRUE)
    targets <- c(targets, merged_dirs)
    labels <- c(labels, paste("Merged:", basename(merged_dirs)))
  }
  
  stats::setNames(targets, labels)
}

list_files_by_pattern <- function(root_dir, pattern) {
  if (is.null(root_dir) || !dir.exists(root_dir)) {
    return(character())
  }
  list.files(root_dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}

build_cli_args <- function(input) {
  args <- c("--config", normalizePath(input$config_path, winslash = "/", mustWork = FALSE))
  
  if (nzchar(trimws(input$sample_sheet))) {
    args <- c(args, "--sample-sheet", normalizePath(input$sample_sheet, winslash = "/", mustWork = FALSE))
  }
  if (nzchar(input$mode)) {
    args <- c(args, "--mode", input$mode)
  }
  if (length(input$samples) > 0) {
    args <- c(args, "--samples", paste(input$samples, collapse = ","))
  }
  if (length(input$pairs) > 0) {
    args <- c(args, "--pairs", paste(input$pairs, collapse = ","))
  }
  if (length(input$groups) > 0) {
    args <- c(args, "--groups", paste(input$groups, collapse = ","))
  }
  if (length(input$sample_steps) > 0 && input$mode %in% c("separate", "paired", "all")) {
    args <- c(args, "--sample-steps", paste(input$sample_steps, collapse = ","))
  }
  if (length(input$merged_steps) > 0 && input$mode %in% c("merged", "all")) {
    args <- c(args, "--merged-steps", paste(input$merged_steps, collapse = ","))
  }
  if (nzchar(input$from_step) && input$mode %in% c("separate", "paired", "all")) {
    args <- c(args, "--from-step", input$from_step)
  }
  if (nzchar(input$merged_from_step) && input$mode %in% c("merged", "all")) {
    args <- c(args, "--merged-from-step", input$merged_from_step)
  }
  if (nzchar(trimws(input$run_label)) && input$mode %in% c("merged", "all")) {
    args <- c(args, "--run-label", trimws(input$run_label))
  }
  
  args
}

rscript_binary <- function() {
  file.path(R.home("bin"), "Rscript")
}

ui <- page_sidebar(
  title = div(class = "portal-title", "CosMx Kidney Portal"),
  theme = bs_theme(
    version = 5,
    bg = "#f4ede3",
    fg = "#1f1b17",
    primary = "#8a3643",
    secondary = "#1f5b59"
  ),
  sidebar = sidebar(
    width = 360,
    textInput("config_path", "Config file", value = file.path(repo_dir, "config.yaml")),
    textInput("sample_sheet", "Sample sheet (optional)", value = ""),
    actionButton("refresh_context", "Refresh samples"),
    selectInput(
      "mode",
      "Run mode",
      choices = c("All" = "all", "Separate" = "separate", "Paired" = "paired", "Merged" = "merged"),
      selected = "all"
    ),
    selectizeInput("samples", "Samples", choices = NULL, multiple = TRUE),
    selectizeInput("pairs", "Pairs", choices = NULL, multiple = TRUE),
    selectizeInput("groups", "Groups", choices = NULL, multiple = TRUE),
    checkboxGroupInput(
      "sample_steps",
      "Sample steps",
      choices = stats::setNames(SAMPLE_STEP_ORDER, vapply(SAMPLE_STEP_ORDER, pretty_step_label, character(1))),
      selected = SAMPLE_STEP_ORDER
    ),
    selectInput(
      "from_step",
      "Resume sample from",
      choices = c("None" = "", stats::setNames(SAMPLE_STEP_ORDER, vapply(SAMPLE_STEP_ORDER, pretty_step_label, character(1)))),
      selected = ""
    ),
    checkboxGroupInput(
      "merged_steps",
      "Merged steps",
      choices = stats::setNames(MERGED_STEP_ORDER, vapply(MERGED_STEP_ORDER, pretty_step_label, character(1))),
      selected = MERGED_STEP_ORDER
    ),
    selectInput(
      "merged_from_step",
      "Resume merged from",
      choices = c("None" = "", stats::setNames(MERGED_STEP_ORDER, vapply(MERGED_STEP_ORDER, pretty_step_label, character(1)))),
      selected = ""
    ),
    textInput("run_label", "Merged label (optional)", value = ""),
    layout_column_wrap(
      widths = c(1/2, 1/2),
      actionButton("preview_run", "Preview"),
      actionButton("run_pipeline", "Run pipeline", class = "btn-primary")
    )
  ),
  tags$style(HTML("
    .portal-title { font-family: 'IBM Plex Sans', 'Aptos', 'Segoe UI', sans-serif; font-weight: 700; letter-spacing: 0.02em; }
    .bslib-sidebar-layout > .sidebar { border-right: 1px solid rgba(31, 27, 23, 0.1); }
    .card { border: 1px solid rgba(31, 27, 23, 0.08); box-shadow: 0 12px 24px rgba(31, 27, 23, 0.06); }
    .metric-value { font-size: 2rem; font-weight: 700; color: #8a3643; }
    .metric-label { text-transform: uppercase; letter-spacing: 0.08em; font-size: 0.75rem; color: #5b544c; }
    .log-box pre { min-height: 360px; background: #1f1b17; color: #f4ede3; padding: 1rem; border-radius: 0.75rem; }
  ")),
  navset_card_tab(
    nav_panel(
      "Dashboard",
      layout_column_wrap(
        width = 1/3,
        card(card_body(div(class = "metric-label", "Valid samples"), div(class = "metric-value", textOutput("metric_samples")))),
        card(card_body(div(class = "metric-label", "Pair groups"), div(class = "metric-value", textOutput("metric_pairs")))),
        card(card_body(div(class = "metric-label", "Latest run"), div(class = "metric-value", textOutput("metric_latest_run"))))
      ),
      card(
        full_screen = TRUE,
        card_header("Selected samples"),
        DTOutput("sample_table")
      ),
      card(
        full_screen = TRUE,
        card_header("Latest pipeline results"),
        DTOutput("latest_results")
      )
    ),
    nav_panel(
      "Runner",
      card(
        card_header("Command preview"),
        verbatimTextOutput("command_preview")
      ),
      card(
        card_header("Run status"),
        textOutput("run_status")
      ),
      card(
        class = "log-box",
        card_header("Live log"),
        verbatimTextOutput("run_log")
      )
    ),
    nav_panel(
      "Results",
      layout_sidebar(
        sidebar = sidebar(
          selectInput("result_target", "Result folder", choices = NULL),
          selectInput("plot_file", "Plot", choices = NULL),
          selectInput("table_file", "Table", choices = NULL)
        ),
        card(card_header("Plot preview"), imageOutput("plot_preview")),
        card(card_header("Table preview"), DTOutput("table_preview"))
      )
    ),
    nav_panel(
      "B Cells",
      layout_sidebar(
        sidebar = sidebar(
          selectInput("bcell_file", "B-cell table", choices = NULL),
          selectInput("bcell_plot", "B-cell plot", choices = NULL)
        ),
        card(card_header("B-cell table"), DTOutput("bcell_table")),
        card(card_header("B-cell plot"), imageOutput("bcell_plot_preview"))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    process = NULL,
    log_file = NULL,
    exit_status = NULL
  )
  
  refresh_context <- reactiveVal(Sys.time())
  observeEvent(input$refresh_context, refresh_context(Sys.time()))
  
  app_context <- reactive({
    refresh_context()
    cfg <- load_config(input$config_path)
    sample_sheet <- if (nzchar(trimws(input$sample_sheet))) input$sample_sheet else NULL
    sample_tbl <- load_sample_table(cfg, sample_sheet_override = sample_sheet)
    list(config = cfg, samples = sample_tbl)
  })
  
  observe({
    ctx <- app_context()
    sample_tbl <- ctx$samples
    
    updateSelectizeInput(session, "samples", choices = sample_tbl$sample_id, server = TRUE)
    updateSelectizeInput(session, "pairs", choices = sort(unique(na.omit(sample_tbl$pair_id))), server = TRUE)
    updateSelectizeInput(session, "groups", choices = sort(unique(na.omit(sample_tbl$group_id))), server = TRUE)
    
    targets <- list_target_directories(ctx$config$paths$output_dir)
    updateSelectInput(session, "result_target", choices = targets)
    
    bcell_files <- list_files_by_pattern(ctx$config$paths$output_dir, "bcell_interactions\\.csv$|cell_proximity_enrichment\\.csv$|celltype_abundance\\.csv$")
    updateSelectInput(session, "bcell_file", choices = stats::setNames(bcell_files, basename(bcell_files)))
  })
  
  selected_samples <- reactive({
    ctx <- app_context()
    select_samples(
      samples = ctx$samples,
      sample_ids = input$samples %||% NULL,
      pair_ids = input$pairs %||% NULL,
      group_ids = input$groups %||% NULL
    )
  })
  
  output$metric_samples <- renderText({
    nrow(selected_samples())
  })
  
  output$metric_pairs <- renderText({
    sample_tbl <- selected_samples()
    n_distinct(na.omit(sample_tbl$pair_id))
  })
  
  output$metric_latest_run <- renderText({
    latest <- list_latest_run_files(app_context()$config$paths$output_dir)$latest_run
    if (is.null(latest)) {
      "none"
    } else {
      basename(latest)
    }
  })
  
  output$sample_table <- renderDT({
    datatable(selected_samples(), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$latest_results <- renderDT({
    latest <- list_latest_run_files(app_context()$config$paths$output_dir)
    if (!is.null(latest$sample_results) && file.exists(latest$sample_results)) {
      return(datatable(readr::read_csv(latest$sample_results, show_col_types = FALSE), options = list(scrollX = TRUE)))
    }
    if (!is.null(latest$merged_results) && file.exists(latest$merged_results)) {
      return(datatable(readr::read_csv(latest$merged_results, show_col_types = FALSE), options = list(scrollX = TRUE)))
    }
    datatable(tibble(message = "No pipeline results found yet."))
  })
  
  cli_args <- reactive({
    build_cli_args(input)
  })
  
  output$command_preview <- renderText({
    paste(c(shQuote(rscript_binary()), shQuote(file.path(repo_dir, "CosMx_pipeline.R")), shQuote(cli_args())), collapse = " ")
  })
  
  observeEvent(input$run_pipeline, {
    ctx <- app_context()
    log_dir <- ensure_dir(file.path(ctx$config$paths$output_dir, "pipeline_runs", "live_logs"))
    rv$log_file <- file.path(log_dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", input$mode, ".log"))
    rv$exit_status <- NULL
    
    rv$process <- processx::process$new(
      command = rscript_binary(),
      args = c(file.path(repo_dir, "CosMx_pipeline.R"), cli_args()),
      stdout = rv$log_file,
      stderr = rv$log_file,
      supervise = TRUE,
      cleanup_tree = TRUE
    )
  })
  
  observeEvent(input$preview_run, {
    showModal(modalDialog(
      title = "Pipeline preview",
      tags$p(paste("Mode:", input$mode)),
      tags$p(paste("Selected samples:", nrow(selected_samples()))),
      tags$pre(paste(c(rscript_binary(), shQuote(file.path(repo_dir, "CosMx_pipeline.R")), shQuote(cli_args())), collapse = " ")),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observe({
    invalidateLater(1500, session)
    if (!is.null(rv$process) && !rv$process$is_alive() && is.null(rv$exit_status)) {
      rv$exit_status <- rv$process$get_exit_status()
      refresh_context(Sys.time())
    }
  })
  
  output$run_status <- renderText({
    if (is.null(rv$process)) {
      return("No run has been started in this session.")
    }
    if (rv$process$is_alive()) {
      return("Pipeline is running.")
    }
    paste("Pipeline finished with exit status", rv$exit_status %||% 0)
  })
  
  output$run_log <- renderText({
    tail_text(rv$log_file)
  })
  
  observe({
    target_dir <- input$result_target
    plot_files <- list_files_by_pattern(target_dir, "\\.(png|jpg|jpeg)$")
    table_files <- list_files_by_pattern(target_dir, "\\.(csv|tsv|txt)$")
    
    updateSelectInput(session, "plot_file", choices = stats::setNames(plot_files, basename(plot_files)))
    updateSelectInput(session, "table_file", choices = stats::setNames(table_files, basename(table_files)))
  })
  
  output$plot_preview <- renderImage({
    req(input$plot_file)
    list(src = normalizePath(input$plot_file, winslash = "/", mustWork = TRUE))
  }, deleteFile = FALSE)
  
  output$table_preview <- renderDT({
    req(input$table_file)
    datatable(read_table_preview(input$table_file), options = list(scrollX = TRUE, pageLength = 12))
  })
  
  observe({
    req(input$bcell_file)
    bcell_dir <- dirname(input$bcell_file)
    bcell_plots <- list_files_by_pattern(bcell_dir, "cell_proximity_heatmap\\.png$|cell_proximity_network\\.png$")
    updateSelectInput(session, "bcell_plot", choices = stats::setNames(bcell_plots, basename(bcell_plots)))
  })
  
  output$bcell_table <- renderDT({
    req(input$bcell_file)
    datatable(read_table_preview(input$bcell_file), options = list(scrollX = TRUE, pageLength = 12))
  })
  
  output$bcell_plot_preview <- renderImage({
    req(input$bcell_plot)
    list(src = normalizePath(input$bcell_plot, winslash = "/", mustWork = TRUE))
  }, deleteFile = FALSE)
}

shinyApp(ui, server)
