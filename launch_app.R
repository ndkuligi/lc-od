args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
  app_dir <- dirname(script_path)
} else {
  app_dir <- getwd()
}

setwd(app_dir)
Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "TRUE")

source(file.path(app_dir, "renv", "activate.R"))

options(shiny.launch.browser = TRUE)

port <- httpuv::randomPort()
shiny::runApp(appDir = app_dir, host = "127.0.0.1", port = port)