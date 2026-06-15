# LC-OD - Shiny application for the layer charge determination using the OD method
# Based on:
# Kuligiewicz et al. (2015), https://doi.org/10.1346/ccmn.2015.0630603 
# Copyright © Institute of Geological Sciences, Polish Academy of Sciences
# Author: Artur Kuligiewicz
# Version: 1.1.0
# Repository: https://github.com/ndkuligi/lc-od
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(signal)
library(writexl)
library(plotly)

# =========================================================
# Helper functions
# =========================================================

read_ftir_file <- function(file) {
  lines <- readLines(file, n = 20, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  if (length(lines) == 0) {
    stop(paste("File", basename(file), "is empty."))
  }
  
  is_numeric_string <- function(x) {
    grepl(
      pattern = "^[-+]?(?:\\d+(?:[\\.,]\\d*)?|[\\.,]\\d+)(?:[eE][-+]?\\d+)?$",
      x = x
    )
  }
  
  split_line <- function(line, sep) {
    if (sep == "space") {
      parts <- strsplit(line, "\\s+")[[1]]
    } else {
      parts <- strsplit(line, split = sep, fixed = TRUE)[[1]]
    }
    parts <- trimws(parts)
    parts[parts != ""]
  }
  
  candidate_separators <- c(";", ",", "\t", "space")
  
  separator_score <- sapply(candidate_separators, function(sep) {
    score <- 0
    
    for (line in lines) {
      parts <- split_line(line, sep)
      if (length(parts) == 2 && all(is_numeric_string(parts))) {
        score <- score + 1
      }
    }
    
    score
  })
  
  best_sep <- candidate_separators[which.max(separator_score)]
  
  if (max(separator_score) == 0) {
    stop(paste("Could not detect the column separator in file:", basename(file)))
  }
  
  if (best_sep == "space") {
    data <- read.table(
      file = file,
      header = FALSE,
      stringsAsFactors = FALSE,
      strip.white = TRUE,
      fill = TRUE,
      colClasses = "character"
    )
  } else {
    data <- read.table(
      file = file,
      sep = best_sep,
      header = FALSE,
      stringsAsFactors = FALSE,
      strip.white = TRUE,
      colClasses = "character"
    )
  }
  
  if (ncol(data) < 2) {
    stop(paste("File", basename(file), "does not contain at least 2 columns."))
  }
  
  data <- data[, 1:2, drop = FALSE]
  
  data[] <- lapply(data, function(x) {
    x <- trimws(x)
    x <- gsub(",", ".", x, fixed = TRUE)
    as.numeric(x)
  })
  
  colnames(data) <- c("wavenumber", "signal")
  
  data <- data[complete.cases(data), , drop = FALSE]
  
  if (nrow(data) == 0) {
    stop(paste("File", basename(file), "is empty or does not contain valid data."))
  }
  
  data <- data[order(data$wavenumber, decreasing = TRUE), , drop = FALSE]
  
  return(data)
}

prepare_data <- function(files,
                         file_names = NULL,
                         wn_plot_min = 2550,
                         wn_plot_max = 2850,
                         wn_min = 2600,
                         wn_max = 2720,
                         sg_window = 13) {
  
  if (length(files) == 0) {
    stop("No files were provided.")
  }
  
  if (is.null(file_names)) {
    file_names <- basename(files)
  }
  
  # ---------------------------------------------------------
  # Function cleaning spectra
  # ---------------------------------------------------------
  clean_spectrum <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df <- df[, c("wavenumber", "signal"), drop = FALSE]
    df <- df[is.finite(df$wavenumber) & is.finite(df$signal), , drop = FALSE]
    
    if (nrow(df) < 3) return(NULL)
    
    # approx() requires ascending X values
    df <- df[order(df$wavenumber, decreasing = FALSE), , drop = FALSE]
    
    # remove duplicated wavenumber values
    df <- df[!duplicated(df$wavenumber), , drop = FALSE]
    
    if (nrow(df) < 3) return(NULL)
    
    return(df)
  }
  
  # ---------------------------------------------------------
  # Reading files, checking for errors
  # ---------------------------------------------------------
  raw_data_list <- vector("list", length(files))
  skipped_messages <- character(0)
  
  for (i in seq_along(files)) {
    raw_data_list[[i]] <- tryCatch(
      read_ftir_file(files[i]),
      error = function(e) {
        skipped_messages <<- c(
          skipped_messages,
          paste0(file_names[i], " — skipped: ", conditionMessage(e))
        )
        NULL
      }
    )
  }
  
  # ---------------------------------------------------------
  # Find the first correct reference spectrum
  # ---------------------------------------------------------
  reference_index <- NA_integer_
  reference_df <- NULL
  
  for (i in seq_along(raw_data_list)) {
    if (!is.null(raw_data_list[[i]])) {
      cleaned <- clean_spectrum(raw_data_list[[i]])
      if (!is.null(cleaned)) {
        reference_index <- i
        reference_df <- cleaned
        break
      } else {
        skipped_messages <- c(
          skipped_messages,
          paste0(file_names[i], " — skipped: not enough valid points after cleaning.")
        )
      }
    }
  }
  
  if (is.na(reference_index) || is.null(reference_df)) {
    stop("No valid reference file is available after reading and cleaning the data.")
  }
  
  reference_wn <- reference_df$wavenumber
  
  if (length(reference_wn) < 3 || anyNA(reference_wn) || any(!is.finite(reference_wn))) {
    stop("Reference file has invalid wavenumber values.")
  }
  
  # ---------------------------------------------------------
  # Interpolation of all spectra to a common wavelength grid
  # ---------------------------------------------------------
  aligned_signals <- list()
  valid_names <- character(0)
  valid_paths <- character(0)
  
  ref_range <- range(reference_wn, na.rm = TRUE, finite = TRUE)
  
  for (i in seq_along(raw_data_list)) {
    current_raw <- raw_data_list[[i]]
    
    if (is.null(current_raw)) {
      next
    }
    
    current <- clean_spectrum(current_raw)
    
    if (is.null(current)) {
      skipped_messages <- c(
        skipped_messages,
        paste0(file_names[i], " — skipped: not enough valid points after cleaning.")
      )
      next
    }
    
    current_wn <- current$wavenumber
    current_signal <- current$signal
    
    if (length(current_wn) < 3 || anyNA(current_wn) || any(!is.finite(current_wn))) {
      skipped_messages <- c(
        skipped_messages,
        paste0(file_names[i], " — skipped: invalid wavenumber values.")
      )
      next
    }
    
    if (length(current_signal) < 3 || anyNA(current_signal) || any(!is.finite(current_signal))) {
      skipped_messages <- c(
        skipped_messages,
        paste0(file_names[i], " — skipped: invalid signal values.")
      )
      next
    }
    
    cur_range <- range(current_wn, na.rm = TRUE, finite = TRUE)
    
    has_valid_range <- length(cur_range) == 2 &&
      !anyNA(cur_range) &&
      all(is.finite(cur_range))
    
    if (!isTRUE(has_valid_range)) {
      skipped_messages <- c(
        skipped_messages,
        paste0(file_names[i], " — skipped: invalid wavenumber range.")
      )
      next
    }
    
    covers_reference <- (cur_range[1] <= ref_range[1]) && (cur_range[2] >= ref_range[2])
    
    if (!isTRUE(covers_reference)) {
      skipped_messages <- c(
        skipped_messages,
        paste0(file_names[i], " — skipped: does not fully cover the reference wavenumber range.")
      )
      next
    }
    
    same_grid <- isTRUE(all.equal(current_wn, reference_wn, tolerance = 1e-4))
    
    if (same_grid) {
      aligned_y <- current_signal
    } else {
      interp <- tryCatch(
        approx(
          x = current_wn,
          y = current_signal,
          xout = reference_wn,
          method = "linear",
          rule = 1,
          ties = mean
        ),
        error = function(e) NULL
      )
      
      if (is.null(interp) || is.null(interp$y)) {
        skipped_messages <- c(
          skipped_messages,
          paste0(file_names[i], " — skipped: interpolation failed.")
        )
        next
      }
      
      aligned_y <- interp$y
    }
    
    if (length(aligned_y) != length(reference_wn)) {
      skipped_messages <- c(
        skipped_messages,
        paste0(file_names[i], " — skipped: interpolated signal has wrong length.")
      )
      next
    }
    
    if (anyNA(aligned_y) || any(!is.finite(aligned_y))) {
      skipped_messages <- c(
        skipped_messages,
        paste0(file_names[i], " — skipped: interpolated signal contains NA/Inf.")
      )
      next
    }
    
    sample_name <- tools::file_path_sans_ext(file_names[i])
    
    aligned_signals[[sample_name]] <- aligned_y
    valid_names <- c(valid_names, file_names[i])
    valid_paths <- c(valid_paths, files[i])
  }
  
  if (length(aligned_signals) == 0) {
    stop("No valid spectra available after cleaning and interpolation.")
  }
  
  # ---------------------------------------------------------
  # Results table
  # ---------------------------------------------------------
  result <- data.frame(wavenumber = reference_wn)
  
  for (nm in names(aligned_signals)) {
    result[[nm]] <- aligned_signals[[nm]]
  }
  
  # ---------------------------------------------------------
  # X axis step
  # ---------------------------------------------------------
  delta_wn <- abs(mean(diff(result$wavenumber), na.rm = TRUE))
  
  if (!is.finite(delta_wn) || is.na(delta_wn) || delta_wn <= 0) {
    stop("Invalid wavenumber spacing detected.")
  }
  
  # ---------------------------------------------------------
  # Second derivative
  # ---------------------------------------------------------
  result_d2 <- data.frame(wavenumber = result$wavenumber)
  
  for (i in 2:ncol(result)) {
    current_y <- result[[i]]
    
    if (length(current_y) < sg_window || anyNA(current_y) || any(!is.finite(current_y))) {
      result_d2[[colnames(result)[i]]] <- rep(NA_real_, length(current_y))
    } else {
      result_d2[[colnames(result)[i]]] <- tryCatch(
        signal::sgolayfilt(
          x = current_y,
          p = 3,
          n = sg_window,
          m = 2,
          ts = delta_wn
        ),
        error = function(e) rep(NA_real_, length(current_y))
      )
    }
  }
  
  result_d2_range <- dplyr::filter(
    result_d2,
    wavenumber >= wn_plot_min & wavenumber <= wn_plot_max
  )
  
  # ---------------------------------------------------------
  # Finding minima
  # ---------------------------------------------------------
  peak_positions <- data.frame(
    sample = character(),
    peak_wavenumber_raw = numeric(),
    peak_wavenumber_fit = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 2:ncol(result_d2)) {
    sample_name <- colnames(result_d2)[i]
    
    range_data <- result_d2[
      result_d2$wavenumber >= wn_min & result_d2$wavenumber <= wn_max,
      c("wavenumber", sample_name),
      drop = FALSE
    ]
    
    colnames(range_data) <- c("wavenumber", "d2")
    range_data <- range_data[
      is.finite(range_data$wavenumber) & is.finite(range_data$d2),
      ,
      drop = FALSE
    ]
    
    if (nrow(range_data) < 5) {
      peak_raw <- NA_real_
      peak_fit <- NA_real_
    } else {
      idx_min <- which.min(range_data$d2)
      
      if (length(idx_min) == 0 || is.na(idx_min)) {
        peak_raw <- NA_real_
        peak_fit <- NA_real_
      } else {
        peak_raw <- range_data$wavenumber[idx_min]
        
        if (idx_min <= 2 || idx_min >= (nrow(range_data) - 1)) {
          peak_fit <- NA_real_
        } else {
          fragment <- range_data[(idx_min - 2):(idx_min + 2), , drop = FALSE]
          
          if (nrow(fragment) < 5 || anyNA(fragment) || any(!is.finite(as.matrix(fragment)))) {
            peak_fit <- NA_real_
          } else {
            model <- tryCatch(
              lm(d2 ~ wavenumber + I(wavenumber^2), data = fragment),
              error = function(e) NULL
            )
            
            if (is.null(model)) {
              peak_fit <- NA_real_
            } else {
              coef_model <- coef(model)
              a <- unname(coef_model["I(wavenumber^2)"])
              b <- unname(coef_model["wavenumber"])
              
              if (length(a) == 0 || length(b) == 0 ||
                  is.na(a) || is.na(b) ||
                  !is.finite(a) || !is.finite(b) ||
                  a == 0) {
                peak_fit <- NA_real_
              } else {
                peak_fit <- -b / (2 * a)
                if (!is.finite(peak_fit)) {
                  peak_fit <- NA_real_
                }
              }
            }
          }
        }
      }
    }
    
    peak_positions <- rbind(
      peak_positions,
      data.frame(
        sample = sample_name,
        peak_wavenumber_raw = peak_raw,
        peak_wavenumber_fit = peak_fit,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # ---------------------------------------------------------
  # Statistics
  # n = number of correct peak_wavenumber_fit used in calculations
  # ---------------------------------------------------------
  positions <- peak_positions$peak_wavenumber_fit
  positions <- positions[is.finite(positions)]
  
  mean_position <- if (length(positions) > 0) mean(positions, na.rm = TRUE) else NA_real_
  sd_position <- if (length(positions) > 1) sd(positions, na.rm = TRUE) else NA_real_
  n_positions <- sum(is.finite(positions))
  
  LC_SFM <- if (is.finite(mean_position)) {
    0.6 - (mean_position - 2686) * 0.03
  } else {
    NA_real_
  }
  
  LC_AAM <- if (is.finite(mean_position)) {
    0.38 - (mean_position - 2686) * 0.015
  } else {
    NA_real_
  }
  
  if (is.finite(mean_position) && is.finite(sd_position) && n_positions > 1) {
    t_crit <- qt(0.975, n_positions - 1)
    sem_position <- sd_position / sqrt(n_positions)
    
    LC_SFM_ERROR <- sqrt((0.01)^2 + ((mean_position - 2686) * 0.001)^2 + (-0.03 * sem_position * t_crit)^2)
    LC_AAM_ERROR <- sqrt((0.01)^2 + ((mean_position - 2686) * 0.001)^2 + (-0.015 * sem_position * t_crit)^2)
  } else {
    LC_SFM_ERROR <- NA_real_
    LC_AAM_ERROR <- NA_real_
  }
  
  band_statistics <- data.frame(
    mean_wavenumber = mean_position,
    sd_wavenumber = sd_position,
    n = n_positions,
    LC_SFM = LC_SFM,
    LC_SFM_ERROR = LC_SFM_ERROR,
    LC_AAM = LC_AAM,
    LC_AAM_ERROR = LC_AAM_ERROR,
    SG_window_points = sg_window
  )
  
  list(
    result = result,
    result_d2 = result_d2,
    result_d2_range = result_d2_range,
    peak_positions = peak_positions,
    band_statistics = band_statistics,
    files = valid_paths,
    file_names = valid_names,
    skipped_messages = unique(skipped_messages),
    wn_plot_min = wn_plot_min,
    wn_plot_max = wn_plot_max,
    wn_min = wn_min,
    wn_max = wn_max,
    sg_window = sg_window
  )
}

build_hover_text <- function(file_name, wn, y, y_label) {
  paste0(
    "File: ", file_name,
    "<br>Wavenumber: ", sprintf("%.2f", wn), " cm^-1",
    "<br>", y_label, ": ", signif(y, 6)
  )
}

# =========================================================
# UI
# =========================================================

ui <- fluidPage(
  titlePanel("LC-OD: Layer Charge calculator using the OD method"),
  fluidRow(
    column(
      width = 4,
      wellPanel(
        h4("Loaded files"),
        fileInput(
          inputId = "files_csv",
          label = "Load spectral files",
          multiple = TRUE,
          accept = c(".csv", ".xy", ".dpt")
        ),
        actionButton("clear_data", "Load new set / clear data"),
        tags$hr(),
        uiOutput("file_list"),
        tags$hr(),
        h4("Skipped files"),
        uiOutput("skipped_list")
      )
    ),
    column(
      width = 8,
      plotlyOutput("plot_atr", height = "360px")
    )
  ),
  fluidRow(
    column(
      width = 4,
      wellPanel(
        radioButtons(
          inputId = "sg_window",
          label = "Savitzky-Golay window size",
          choices = c("11" = 11, "13" = 13, "15" = 15, "17" = 17),
          selected = 13,
          inline = TRUE
        ),
        tags$hr(),
        h4("Export"),
        downloadButton("export_xlsx", "Export results to .xlsx"),
        tags$hr(),
        h4("Statistics"),
        tableOutput("statistics")
      )
    ),
    column(
      width = 8,
      plotlyOutput("plot_d2", height = "360px")
    )
  )
)

# =========================================================
# Server
# =========================================================

server <- function(input, output, session) {
  app_data <- reactiveVal(NULL)
  selected_sample <- reactiveVal(NULL)
  
  load_and_process <- function() {
    req(input$files_csv)
    
    new_data <- prepare_data(
      files = input$files_csv$datapath,
      file_names = input$files_csv$name,
      wn_plot_min = 2550,
      wn_plot_max = 2850,
      wn_min = 2600,
      wn_max = 2720,
      sg_window = as.numeric(input$sg_window)
    )
    
    app_data(new_data)
    selected_sample(NULL)
    
    showNotification(
      paste("Loaded", length(new_data$file_names), "valid file(s) for analysis."),
      type = "message"
    )
    
    if (length(new_data$skipped_messages) > 0) {
      for (msg in new_data$skipped_messages) {
        showNotification(msg, type = "warning", duration = 8)
      }
    }
  }
  
  observeEvent(input$files_csv, {
    req(input$files_csv)
    
    tryCatch(
      load_and_process(),
      error = function(e) {
        app_data(NULL)
        selected_sample(NULL)
        showNotification(conditionMessage(e), type = "error", duration = NULL)
      }
    )
  })
  
  observeEvent(input$sg_window, {
    req(input$files_csv)
    
    tryCatch(
      load_and_process(),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
      }
    )
  }, ignoreInit = TRUE)
  
  observeEvent(input$clear_data, {
    app_data(NULL)
    selected_sample(NULL)
    showNotification(
      "Current data have been cleared. You can now load a new set of files.",
      type = "message"
    )
  })
  
  observeEvent(event_data("plotly_click", source = "atr_plot"), {
    click_data <- event_data("plotly_click", source = "atr_plot")
    
    if (!is.null(click_data) && !is.null(click_data$curveNumber)) {
      data <- app_data()
      req(data)
      
      samples <- colnames(data$result)[-1]
      idx <- click_data$curveNumber + 1
      
      if (idx >= 1 && idx <= length(samples)) {
        clicked_sample <- samples[idx]
        
        if (identical(selected_sample(), clicked_sample)) {
          selected_sample(NULL)
        } else {
          selected_sample(clicked_sample)
        }
      }
    }
  })
  
  output$file_list <- renderUI({
    data <- app_data()
    
    if (is.null(data)) {
      return(tags$p("No files loaded."))
    }
    
    tags$div(
      tags$p(strong("Number of valid files:"), length(data$file_names)),
      tags$ul(lapply(data$file_names, tags$li))
    )
  })
  
  output$skipped_list <- renderUI({
    data <- app_data()
    
    if (is.null(data) || length(data$skipped_messages) == 0) {
      return(tags$p("No skipped files."))
    }
    
    tags$ul(
      lapply(data$skipped_messages, function(msg) {
        tags$li(style = "color:#b22222;", msg)
      })
    )
  })
  
  output$plot_atr <- renderPlotly({
    data <- app_data()
    req(data)
    
    samples <- colnames(data$result)[-1]
    selected <- selected_sample()
    
    if (length(samples) == 0) {
      return(
        plot_ly(
          x = numeric(0),
          y = numeric(0),
          type = "scatter",
          mode = "lines",
          source = "atr_plot"
        ) %>%
          layout(
            xaxis = list(title = "Wavenumber (cm^-1)", autorange = "reversed"),
            yaxis = list(title = "Absorbance"),
            annotations = list(
              text = "No valid spectra to display.",
              x = 0.5, y = 0.5,
              xref = "paper", yref = "paper",
              showarrow = FALSE
            )
          ) %>%
          config(scrollZoom = TRUE, displaylogo = FALSE)
      )
    }
    
    first_sample <- samples[1]
    first_width <- if (!is.null(selected) && identical(first_sample, selected)) 3 else 1.2
    
    p <- plot_ly(
      x = data$result$wavenumber,
      y = data$result[[first_sample]],
      type = "scatter",
      mode = "lines",
      name = first_sample,
      text = build_hover_text(
        file_name = first_sample,
        wn = data$result$wavenumber,
        y = data$result[[first_sample]],
        y_label = "Absorbance"
      ),
      hovertemplate = "%{text}<extra></extra>",
      line = list(width = first_width),
      opacity = 0.85,
      showlegend = FALSE,
      source = "atr_plot"
    )
    
    if (length(samples) > 1) {
      for (sample in samples[-1]) {
        line_width <- if (!is.null(selected) && identical(sample, selected)) 3 else 1.2
        
        p <- add_lines(
          p,
          x = data$result$wavenumber,
          y = data$result[[sample]],
          name = sample,
          text = build_hover_text(
            file_name = sample,
            wn = data$result$wavenumber,
            y = data$result[[sample]],
            y_label = "Absorbance"
          ),
          hovertemplate = "%{text}<extra></extra>",
          line = list(width = line_width),
          opacity = 0.85,
          showlegend = FALSE
        )
      }
    }
    
    p %>%
      layout(
        xaxis = list(
          title = "Wavenumber (cm^-1)",
          autorange = "reversed"
        ),
        yaxis = list(title = "Absorbance"),
        dragmode = "pan"
      ) %>%
      config(scrollZoom = TRUE, displaylogo = FALSE)
  })
  
  output$plot_d2 <- renderPlotly({
    data <- app_data()
    req(data)
    
    samples <- colnames(data$result_d2_range)[-1]
    selected <- selected_sample()
    minima_df <- data$peak_positions %>%
      dplyr::filter(!is.na(peak_wavenumber_fit), is.finite(peak_wavenumber_fit))
    
    if (length(samples) == 0) {
      return(
        plot_ly(
          x = numeric(0),
          y = numeric(0),
          type = "scatter",
          mode = "lines",
          source = "d2_plot"
        ) %>%
          layout(
            xaxis = list(
              title = "Wavenumber (cm^-1)",
              autorange = "reversed",
              range = c(data$wn_plot_max, data$wn_plot_min)
            ),
            yaxis = list(title = "Second derivative"),
            annotations = list(
              text = "No valid second-derivative data to display.",
              x = 0.5, y = 0.5,
              xref = "paper", yref = "paper",
              showarrow = FALSE
            )
          ) %>%
          config(scrollZoom = TRUE, displaylogo = FALSE)
      )
    }
    
    y_values <- as.matrix(data$result_d2_range[, -1, drop = FALSE])
    y_values <- y_values[is.finite(y_values)]
    
    if (length(y_values) == 0) {
      y_min <- -1
      y_max <- 1
    } else {
      y_min <- min(y_values, na.rm = TRUE)
      y_max <- max(y_values, na.rm = TRUE)
    }
    
    first_sample <- samples[1]
    first_width <- if (!is.null(selected) && identical(first_sample, selected)) 3 else 1.2
    
    p <- plot_ly(
      x = data$result_d2_range$wavenumber,
      y = data$result_d2_range[[first_sample]],
      type = "scatter",
      mode = "lines",
      name = first_sample,
      text = build_hover_text(
        file_name = first_sample,
        wn = data$result_d2_range$wavenumber,
        y = data$result_d2_range[[first_sample]],
        y_label = "Second derivative"
      ),
      hovertemplate = "%{text}<extra></extra>",
      line = list(width = first_width),
      opacity = 0.85,
      showlegend = FALSE,
      source = "d2_plot"
    )
    
    if (length(samples) > 1) {
      for (sample in samples[-1]) {
        line_width <- if (!is.null(selected) && identical(sample, selected)) 3 else 1.2
        
        p <- add_lines(
          p,
          x = data$result_d2_range$wavenumber,
          y = data$result_d2_range[[sample]],
          name = sample,
          text = build_hover_text(
            file_name = sample,
            wn = data$result_d2_range$wavenumber,
            y = data$result_d2_range[[sample]],
            y_label = "Second derivative"
          ),
          hovertemplate = "%{text}<extra></extra>",
          line = list(width = line_width),
          opacity = 0.85,
          showlegend = FALSE
        )
      }
    }
    
    if (nrow(minima_df) > 0) {
      for (i in seq_len(nrow(minima_df))) {
        sample_name <- minima_df$sample[i]
        peak_x <- minima_df$peak_wavenumber_fit[i]
        
        if (!is.finite(peak_x)) next
        
        selected_width <- if (!is.null(selected) && identical(sample_name, selected)) 3 else 1.2
        
        p <- add_segments(
          p,
          x = peak_x,
          xend = peak_x,
          y = y_min,
          yend = y_max,
          inherit = FALSE,
          line = list(width = selected_width, dash = "dash"),
          hoverinfo = "text",
          text = paste0(
            "File: ", sample_name,
            "<br>Peak position: ", sprintf("%.2f", peak_x), " cm^-1"
          ),
          showlegend = FALSE
        )
      }
    }
    
    p %>%
      layout(
        xaxis = list(
          title = "Wavenumber (cm^-1)",
          autorange = "reversed",
          range = c(data$wn_plot_max, data$wn_plot_min)
        ),
        yaxis = list(title = "Second derivative"),
        dragmode = "pan"
      ) %>%
      config(scrollZoom = TRUE, displaylogo = FALSE)
  })
  
  output$statistics <- renderTable({
    data <- app_data()
    req(data)
    
    stats_display <- data$band_statistics %>%
      dplyr::select(-SG_window_points, -LC_SFM_ERROR, -LC_AAM_ERROR)
    
    stats_display
  }, digits = 4)
  
  output$export_xlsx <- downloadHandler(
    filename = function() {
      data <- app_data()
      req(data)
      
      first_name <- if (length(data$file_names) > 0) {
        tools::file_path_sans_ext(basename(data$file_names[1]))
      } else {
        "results"
      }
      
      paste0(first_name, "_ROD.xlsx")
    },
    content = function(file) {
      data <- app_data()
      req(data)
      
      skipped_df <- if (length(data$skipped_messages) > 0) {
        data.frame(message = data$skipped_messages, stringsAsFactors = FALSE)
      } else {
        data.frame(message = "No skipped files.", stringsAsFactors = FALSE)
      }
      
      write_xlsx(
        list(
          spectra = data$result,
          second_derivative = data$result_d2,
          peak_positions = data$peak_positions,
          statistics = data$band_statistics,
          skipped_files = skipped_df
        ),
        path = file
      )
    }
  )
}

shinyApp(ui = ui, server = server)
