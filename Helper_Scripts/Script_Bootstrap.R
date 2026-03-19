bootstrap_pipeline_environment <- function(repo_dir,
                                           load_pipeline_utils = FALSE,
                                           verbose = TRUE) {
  repo_dir <- normalizePath(repo_dir, winslash = "/", mustWork = FALSE)
  helper_dir <- file.path(repo_dir, "Helper_Scripts")
  
  Sys.setenv(
    COSMX_PROJECT_DIR = repo_dir,
    COSMX_HELPER_DIR = helper_dir
  )
  
  setup_script <- file.path(repo_dir, "00_Setup.R")
  if (!exists(".cosmx_env_setup_complete", envir = .GlobalEnv)) {
    if (!file.exists(setup_script)) {
      stop("Setup script not found: ", setup_script)
    }
    source(setup_script, local = .GlobalEnv)
  } else if (isTRUE(verbose)) {
    message("CosMx environment already initialized")
  }
  
  if (isTRUE(load_pipeline_utils)) {
    pipeline_utils <- file.path(helper_dir, "Pipeline_Utils.R")
    if (file.exists(pipeline_utils) &&
        !exists("save_giotto_checkpoint", envir = .GlobalEnv)) {
      source(pipeline_utils, local = .GlobalEnv)
    }
  }
  
  invisible(
    list(
      repo_dir = repo_dir,
      helper_dir = helper_dir
    )
  )
}
