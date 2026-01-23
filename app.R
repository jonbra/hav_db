# =============================================================================
# HAV Database - Shiny Application
# =============================================================================
# A web interface for managing sequences, searching, and running analyses.
library(DT)
library(shiny)
library(shinydashboard)
library(plotly)
library(DBI)
library(RSQLite)
library(Biostrings)
## 'msa' package is not required; external 'mafft' is used instead.
library(ape)
library(ggplot2)
library(httr)
library(jsonlite)

# Source helper functions
`%||%` <- function(x, y) if (is.null(x)) y else x

# Determine the application directory in a robust way. Accessing
# `sys.frame(1)$ofile` can fail when running non-interactively (error:
# "not that many frames on the stack"). Try multiple strategies and fall
# back to the current working directory.
app_dir <- (function() {
  # 1) Try to read from sys.frame(1)$ofile if available
  dir <- NULL
  try({
    f <- sys.frame(1)
    if (!is.null(f$ofile)) dir <- dirname(f$ofile)
  }, silent = TRUE)

  # 2) Try commandArgs to detect --file= when invoked via Rscript
  if (is.null(dir) || dir == "") {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      path <- sub("^--file=", "", file_arg[1])
      dir <- dirname(path)
    }
  }

  # 3) Fallback to getwd()
  if (is.null(dir) || dir == "") dir <- getwd()
  dir
})()

# Load optional config (sets `DB_PATH`, can be overridden via HAV_DB_PATH)
local_cfg <- file.path(app_dir, "config_local.R")
if (file.exists(local_cfg)) {
  source(local_cfg)
} else if (file.exists(file.path(app_dir, "config.R"))) {
  local_cfg <- file.path(app_dir, "config_local.R")
  if (file.exists(local_cfg)) {
    source(local_cfg)
  } else if (file.exists(file.path(app_dir, "config.R"))) {
    source(file.path(app_dir, "config.R"))
  }
}

# Source helper functions
source(file.path(app_dir, "R", "db_functions.R"))
source(file.path(app_dir, "R", "sequence_functions.R"))
source(file.path(app_dir, "R", "analysis_functions.R"))
source(file.path(app_dir, "R", "provenance_helpers.R"))

# Source UI components and build UI
source(file.path(app_dir, "R", "ui_components.R"))
ui <- build_ui(app_dir)

# Source server modules
source(file.path(app_dir, "R", "modules", "analysis_module.R"))
source(file.path(app_dir, "R", "modules", "microreact_module.R"))

# Define server function (the remainder of this file implements the server body)
server <- function(input, output, session) {
  rv <- reactiveValues(
    selected_for_analysis = character(0),
    msa_result = NULL,
    tree_result = NULL,
    microreact_token = NULL,
    microreact_team_id = NULL,
    microreact_projects = data.frame(name = character(0), url = character(0), id = character(0), created = character(0), stringsAsFactors = FALSE),
    current_microreact_url = NULL,
    auto_open_local = FALSE
  )

  con <- reactive({
    if (!exists("DB_PATH") || is.null(DB_PATH) || !file.exists(DB_PATH)) {
      stop("Database not found. Please set DB_PATH in config.R or create the database using setup_database.R")
    }
    DBI::dbConnect(RSQLite::SQLite(), DB_PATH)
  })

  session$onSessionEnded(function() {
    try({
      c <- con()
      if (!is.null(c) && DBI::dbIsValid(c)) DBI::dbDisconnect(c)
    }, silent = TRUE)
  })

  # ==========================
  # Dashboard outputs
  # ==========================
  output$total_sequences <- renderInfoBox({
    req(con())
    stats <- tryCatch(get_db_stats(con()), error = function(e) NULL)
    n <- if (!is.null(stats)) stats$total_sequences else NA
    infoBox("Total sequences", n, icon = icon("dna"), color = "purple")
  })

  output$total_viruses <- renderInfoBox({
    req(con())
    n <- tryCatch(DBI::dbGetQuery(con(), "SELECT COUNT(DISTINCT genotype) as n FROM metadata")$n, error = function(e) NA)
    infoBox("Genotypes", n, icon = icon("th"), color = "green")
  })

  output$date_range <- renderInfoBox({
    req(con())
    dr <- tryCatch(DBI::dbGetQuery(con(), "SELECT MIN(sampling_date) as min_date, MAX(sampling_date) as max_date FROM metadata"), error = function(e) NULL)
    txt <- if (!is.null(dr) && nrow(dr) > 0 && !is.na(dr$min_date)) paste0(dr$min_date, " - ", dr$max_date) else "N/A"
    infoBox("Sampling range", txt, icon = icon("calendar-alt"), color = "blue")
  })

  output$virus_plot <- renderPlotly({
    req(con())
    df <- tryCatch(DBI::dbGetQuery(con(), "SELECT genotype, COUNT(*) AS n FROM metadata GROUP BY genotype"), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    plot_ly(df, x = ~genotype, y = ~n, type = 'bar') %>% layout(xaxis = list(title = 'Genotype'), yaxis = list(title = 'Count'))
  })

  output$location_plot <- renderPlotly({
    req(con())
    df <- tryCatch(DBI::dbGetQuery(con(), "SELECT geo_country, COUNT(*) AS n FROM metadata GROUP BY geo_country"), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    plot_ly(df, x = ~geo_country, y = ~n, type = 'bar') %>% layout(xaxis = list(title = 'Country'), yaxis = list(title = 'Count'))
  })

  output$recent_sequences <- DT::renderDataTable({
    req(con())
    df <- tryCatch(DBI::dbGetQuery(con(), "SELECT s.sample_id, m.sample_year, m.genotype, m.geo_country, s.created_at FROM sequences s LEFT JOIN metadata m ON s.sample_id = m.sample_id ORDER BY s.created_at DESC LIMIT 15"), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    datatable(df, options = list(pageLength = 15), rownames = FALSE)
  })

  # Register server modules
  tryCatch({
    register_analysis_server(input, output, session, rv, con)
  }, error = function(e) {
    showNotification(paste("Error registering analysis module:", e$message), type = "error")
  })

  tryCatch({
    register_microreact_server(input, output, session, rv, con, app_dir)
  }, error = function(e) {
    showNotification(paste("Error registering microreact module:", e$message), type = "error")
  })

  # ============================================
  # Browse
  # ============================================
  
  output$all_sequences_table <- renderDT({
    req(con())
    data <- get_sequences_with_metadata(con())
    # Don't show full sequence in table
    data$sequence <- substr(data$sequence, 1, 50)
    data$sequence <- paste0(data$sequence, "...")
    datatable(data, selection = "single", options = list(pageLength = 10), rownames = FALSE)
  })
  
  output$selected_sequence_info <- renderPrint({
    req(input$all_sequences_table_rows_selected)
    req(con())
    data <- get_sequences_with_metadata(con())
    row <- data[input$all_sequences_table_rows_selected, ]
    cat("Sample ID:", row$sample_id, "\n")
    cat("Length:", row$sequence_length, "bp\n")
    cat("Genotype:", row$genotype %||% "N/A", "\n")
    cat("Variant:", row$variant %||% "N/A", "\n")
    cat("Patient ID:", row$patient_id %||% "N/A", "\n")
    cat("Country:", row$geo_country %||% "N/A", "\n")
    cat("Location:", row$geo_location %||% "N/A", "\n")
    cat("Year:", row$sample_year %||% "N/A", "\n")
  })
  
  output$selected_sequence_seq <- renderPrint({
    req(input$all_sequences_table_rows_selected)
    req(con())
    data <- get_sequences_with_metadata(con())
    seq <- data$sequence[input$all_sequences_table_rows_selected]
    # Format sequence with line breaks
    cat(gsub("(.{80})", "\\1\n", seq))
  })
  
  # ============================================
  # Search
  # ============================================
  
  # Update filter dropdowns
  observe({
    req(con())
    genotypes <- c("All" = "", get_unique_values(con(), "genotype"))
    variants <- c("All" = "", get_unique_values(con(), "variant"))
    countries <- c("All" = "", get_unique_values(con(), "geo_country"))
    locations <- c("All" = "", get_unique_values(con(), "geo_location"))
    
    updateSelectInput(session, "search_genotype", choices = genotypes)
    updateSelectInput(session, "search_variant", choices = variants)
    updateSelectInput(session, "search_country", choices = countries)
    updateSelectInput(session, "search_location", choices = locations)
  })
  
  search_data <- eventReactive(input$search_btn, {
    req(con())
    
    # Debug: print search parameters to console
    cat("\n=== Search triggered ===\n")
    cat("Genotype:", input$search_genotype, "\n")
    cat("Variant:", input$search_variant, "\n")
    cat("Patient ID:", input$search_patient, "\n")
    cat("Country:", input$search_country, "\n")
    cat("Location:", input$search_location, "\n")
    cat("Year from:", input$search_year_from, "\n")
    cat("Year to:", input$search_year_to, "\n")
    cat("Min length:", input$search_min_len, "\n")
    cat("Max length:", input$search_max_len, "\n")
    
    result <- search_sequences(con(),
      genotype = if (input$search_genotype == "") NULL else input$search_genotype,
      variant = if (input$search_variant == "") NULL else input$search_variant,
      patient_id = if (input$search_patient == "") NULL else input$search_patient,
      geo_country = if (input$search_country == "") NULL else input$search_country,
      geo_location = if (input$search_location == "") NULL else input$search_location,
      year_from = if (is.na(input$search_year_from)) NULL else input$search_year_from,
      year_to = if (is.na(input$search_year_to)) NULL else input$search_year_to,
      min_length = if (is.na(input$search_min_len)) NULL else input$search_min_len,
      max_length = if (is.na(input$search_max_len)) NULL else input$search_max_len
    )
    # Apply missing value filters
    filters <- input$search_filter_missing
    if (!is.null(filters) && is.data.frame(result) && nrow(result) > 0) {
      if ("sampling_date" %in% filters) {
        result <- result[!is.na(result$sampling_date) & result$sampling_date != "", ]
      }
      if ("geo_location" %in% filters) {
        result <- result[!is.na(result$geo_location) & result$geo_location != "", ]
      }
      if ("geo_country" %in% filters) {
        result <- result[!is.na(result$geo_country) & result$geo_country != "", ]
      }
      if ("genotype" %in% filters) {
        result <- result[!is.na(result$genotype) & result$genotype != "", ]
      }
    }
    cat("Results found:", if (is.data.frame(result)) nrow(result) else 0, "\n")
    result
  })
  
  output$search_results <- renderDT({
    data <- search_data()
    if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) return(NULL)
    data$sequence <- substr(data$sequence, 1, 30)
    data$sequence <- paste0(data$sequence, "...")
    datatable(data, selection = "multiple", 
              options = list(
                pageLength = 50,
                lengthMenu = c(50, 100, 250, 500),
                scrollX = TRUE
              ), 
              rownames = FALSE)
  })
  
  # Selection count display
  output$selection_count <- renderText({
    n_selected <- length(input$search_results_rows_selected)
    data <- search_data()
    n_total <- if (!is.null(data) && is.data.frame(data)) nrow(data) else 0
    paste0(n_selected, " of ", n_total, " sequences selected")
  })
  
  # Select All button - selects all rows in current search results
  observeEvent(input$select_all_results, {
    data <- search_data()
    if (!is.null(data) && is.data.frame(data) && nrow(data) > 0) {
      proxy <- dataTableProxy("search_results")
      selectRows(proxy, 1:nrow(data))
    }
  })
  
  # Deselect All button
  observeEvent(input$deselect_all_results, {
    proxy <- dataTableProxy("search_results")
    selectRows(proxy, NULL)
  })
  
  observeEvent(input$clear_search, {
    updateSelectInput(session, "search_genotype", selected = "")
    updateSelectInput(session, "search_variant", selected = "")
    updateTextInput(session, "search_patient", value = "")
    updateSelectInput(session, "search_country", selected = "")
    updateSelectInput(session, "search_location", selected = "")
    updateNumericInput(session, "search_year_from", value = NA)
    updateNumericInput(session, "search_year_to", value = NA)
    updateNumericInput(session, "search_min_len", value = NA)
    updateNumericInput(session, "search_max_len", value = NA)
  })
  
  observeEvent(input$select_for_analysis, {
    req(input$search_results_rows_selected)
    data <- search_data()
    selected_ids <- data$sample_id[input$search_results_rows_selected]
    rv$selected_for_analysis <- unique(c(rv$selected_for_analysis, selected_ids))
    showNotification(paste("Added", length(selected_ids), "sequences to analysis"), type = "message")
  })
  
  # ============================================
  # Add Data
  # ============================================
  
  observeEvent(input$add_sequence_btn, {
    req(input$new_sample_id, input$new_sequence)
    req(con())
    
    tryCatch({
      # Extract year from date if provided
      sample_year <- if (!is.na(input$new_date)) as.integer(format(input$new_date, "%Y")) else NA
      
      add_sequence_with_metadata(con(),
        sample_id = input$new_sample_id,
        sequence = input$new_sequence,
        sampling_date = as.character(input$new_date),
        sample_year = sample_year,
        genotype = input$new_genotype,
        variant = input$new_variant,
        patient_id = input$new_patient_id,
        geo_location = input$new_location,
        geo_country = input$new_country,
        source = input$new_source,
        transmission_route = input$new_transmission,
        comment = input$new_comment
      )
      showNotification("Sequence added successfully!", type = "message")
      
      # Clear inputs
      updateTextInput(session, "new_sample_id", value = "")
      updateTextAreaInput(session, "new_sequence", value = "")
      updateTextInput(session, "new_genotype", value = "")
      updateTextInput(session, "new_variant", value = "")
      updateTextInput(session, "new_patient_id", value = "")
      updateTextInput(session, "new_location", value = "")
      updateTextInput(session, "new_country", value = "")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  observeEvent(input$import_btn, {
    req(input$fasta_file)
    req(con())
    
    # Read metadata if provided
    metadata_df <- NULL
    if (!is.null(input$metadata_file)) {
      metadata_df <- read.csv(input$metadata_file$datapath, stringsAsFactors = FALSE)
    }
    
    result <- import_fasta_to_db(con(), 
      input$fasta_file$datapath, 
      metadata_df,
      skip_existing = input$skip_existing
    )
    
    output$import_result <- renderPrint({
      cat("Import complete!\n")
      cat("Imported:", result$imported, "\n")
      cat("Skipped:", result$skipped, "\n")
      if (length(result$errors) > 0) {
        cat("Errors:\n")
        cat(paste(result$errors, collapse = "\n"))
      }
    })
  })
  
  # Analysis handlers and UI are provided by `R/modules/analysis_module.R`.
  # The module is registered earlier via `register_analysis_server(input, output, session, rv, con)`.
  
  # ============================================
  # Export
  # ============================================
  
  output$export_fasta <- downloadHandler(
    filename = function() paste0("sequences_", Sys.Date(), ".fasta"),
    content = function(file) {
      req(con())
      export_to_fasta(con(), file, include_metadata = input$export_include_meta)
    }
  )
  
  output$export_selected_fasta <- downloadHandler(
    filename = function() paste0("selected_sequences_", Sys.Date(), ".fasta"),
    content = function(file) {
      req(con())
      req(length(rv$selected_for_analysis) > 0)
      export_to_fasta(con(), file, rv$selected_for_analysis, input$export_include_meta)
    }
  )

  # Download sequences selected in the Search tab (current search results selection)
  output$download_selected_search <- downloadHandler(
    filename = function() paste0("search_selected_sequences_", Sys.Date(), ".fasta"),
    content = function(file) {
      req(con())
      # Use search_data() reactive from server scope; if not available, query DB directly
      data <- tryCatch(search_data(), error = function(e) NULL)
      sel <- input$search_results_rows_selected
      if (is.null(sel) || length(sel) == 0 || is.null(data) || nrow(data) == 0) {
        # Write placeholder FASTA
        writeLines(c(">no_sequences_selected", ""), file)
        return()
      }
      ids <- data$sample_id[sel]
      export_to_fasta(con(), file, ids, include_metadata = FALSE)
    }
  )
  
  output$export_csv <- downloadHandler(
    filename = function() paste0("metadata_", Sys.Date(), ".csv"),
    content = function(file) {
      req(con())
      export_metadata_csv(con(), file)
    }
  )
  
  output$export_selected_csv <- downloadHandler(
    filename = function() paste0("selected_metadata_", Sys.Date(), ".csv"),
    content = function(file) {
      req(con())
      req(length(rv$selected_for_analysis) > 0)
      export_metadata_csv(con(), file, rv$selected_for_analysis)
    }
  )
  
  # Microreact handlers implemented in R/modules/microreact_module.R

  # Analysis persistence UI and handlers are implemented in the analysis module
}

# =============================================================================
# Run App
# =============================================================================

# Read host/port from environment (useful when deploying on a server).
# Defaults: bind to all interfaces and port 3838 for remote testing.
host <- Sys.getenv("SHINY_HOST", "0.0.0.0")
port <- as.integer(Sys.getenv("SHINY_PORT", "3838"))

# When running interactively (e.g. RStudio), keep the usual `shinyApp()` behavior.
# For non-interactive/server runs, call `shiny::runApp()` with explicit host/port.
if (interactive()) {
  shinyApp(ui = ui, server = server)
} else {
  shiny::runApp(list(ui = ui, server = server), host = host, port = port)
}
