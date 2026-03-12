#!/usr/bin/env Rscript

current_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  }
  ofiles <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else frame$ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles) > 0) {
    return(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(getwd(), "Run_CosMx_Portal.R"), winslash = "/", mustWork = FALSE)
}

script_dir <- dirname(current_script_path())
app_dir <- file.path(script_dir, "Shiny_Portal")
args <- commandArgs(trailingOnly = TRUE)
host <- if (length(args) >= 1) args[[1]] else "127.0.0.1"
port <- if (length(args) >= 2) as.integer(args[[2]]) else 3838L

shiny::runApp(appDir = app_dir, host = host, port = port, launch.browser = FALSE)
