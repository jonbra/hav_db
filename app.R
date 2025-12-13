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
library(msa)
library(ape)
library(ggplot2)
library(httr)
library(jsonlite)

# Source helper functions
app_dir <- dirname(sys.frame(1)$ofile %||% ".")
if (app_dir == ".") app_dir <- getwd()

source(file.path(app_dir, "R", "db_functions.R"))
source(file.path(app_dir, "R", "sequence_functions.R"))
source(file.path(app_dir, "R", "analysis_functions.R"))

# Database path
DB_PATH <- file.path(app_dir, "data", "hav.db")

# Null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# UI
# =============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "HAV Database"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Browse Data", tabName = "browse", icon = icon("database")),
      menuItem("Search", tabName = "search", icon = icon("search")),
      menuItem("Add Data", tabName = "add", icon = icon("plus")),
      menuItem("Analysis", tabName = "analysis", icon = icon("dna")),
      menuItem("Microreact", tabName = "microreact", icon = icon("project-diagram")),
      menuItem("Export", tabName = "export", icon = icon("download")),
      menuItem("Settings", tabName = "settings", icon = icon("cog"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("\
        .sequence-display {\
          font-family: 'Courier New', monospace;\
          font-size: 12px;\
          word-wrap: break-word;\
          background-color: #f5f5f5;\
          padding: 10px;\
          border-radius: 4px;\
          max-height: 200px;\
          overflow-y: auto;\
        }\
        .info-box-icon { background-color: rgba(0,0,0,0.1) !important; }\
        /* Make the left sidebar fixed while scrolling */\
        .main-sidebar {\
          position: fixed !important;\
          top: 50px; /* adjust if header height differs */\
          left: 0;\
          height: calc(100vh - 50px);\
          overflow-y: auto;\
          z-index: 1000;\
        }\
        /* Ensure content area leaves space for fixed sidebar */\
        .content-wrapper, .main-footer {\
          margin-left: 230px !important; /* match sidebar width */\
        }\
        /* Keep header fixed and above other elements */
        .main-header {
          position: fixed !important;
          top: 0;
          left: 0;
          right: 0;
          z-index: 1100;
          width: 100%;
        }
        /* Ensure content is pushed below the fixed header */
        .content-wrapper, .main-footer {
          margin-top: 50px !important; /* same as header height */
        }
        /* Adjust sidebar top to sit below fixed header */
        .main-sidebar {
          top: 50px; /* adjust if header height differs */
        }
      "))
    ),
    
    tabItems(
      # Dashboard Tab
      tabItem(tabName = "dashboard",
        fluidRow(
          infoBoxOutput("total_sequences", width = 4),
          infoBoxOutput("total_viruses", width = 4),
          infoBoxOutput("date_range", width = 4)
        ),
        fluidRow(
          box(title = "Sequences by Virus", status = "primary", solidHeader = TRUE,
              plotlyOutput("virus_plot", height = 300), width = 6),
          box(title = "Sequences by Location", status = "primary", solidHeader = TRUE,
              plotlyOutput("location_plot", height = 300), width = 6)
        ),
        fluidRow(
          box(title = "Recent Sequences", status = "info", solidHeader = TRUE,
              DT::dataTableOutput("recent_sequences"), width = 12)
        )
      ),
      
      # Browse Tab
      tabItem(tabName = "browse",
        fluidRow(
          box(title = "All Sequences", status = "primary", solidHeader = TRUE, width = 12,
              DT::dataTableOutput("all_sequences_table"),
              hr(),
              h4("Selected Sequence Details"),
              verbatimTextOutput("selected_sequence_info"),
              div(class = "sequence-display", verbatimTextOutput("selected_sequence_seq"))
          )
        )
      ),
      
      # Search Tab
      tabItem(tabName = "search",
        fluidRow(
          box(title = "Search Filters", status = "warning", solidHeader = TRUE, width = 4,
              selectInput("search_genotype", "Genotype:", choices = c("All" = ""), selected = ""),
              selectInput("search_variant", "Variant:", choices = c("All" = ""), selected = ""),
              textInput("search_patient", "Patient ID:", placeholder = "Search patient ID..."),
              selectInput("search_country", "Country:", choices = c("All" = ""), selected = ""),
              selectInput("search_location", "Location:", choices = c("All" = ""), selected = ""),
              numericInput("search_year_from", "Year From:", value = NA, min = 1990, max = 2100),
              numericInput("search_year_to", "Year To:", value = NA, min = 1990, max = 2100),
              numericInput("search_min_len", "Min Length:", value = NA, min = 0),
              numericInput("search_max_len", "Max Length:", value = NA, min = 0),
              actionButton("search_btn", "Search", class = "btn-primary", width = "100%"),
              hr(),
              actionButton("clear_search", "Clear Filters", width = "100%")
          ),
          checkboxGroupInput("search_filter_missing", "Hide samples missing:",
            choices = list(
              "Sampling date" = "sampling_date",
              "Country" = "geo_country",
              "Genotype" = "genotype"
            ),
            inline = FALSE
          ),
          box(title = "Search Results", status = "primary", solidHeader = TRUE, width = 8,
              DT::dataTableOutput("search_results"),
              hr(),
              fluidRow(
                column(4, actionButton("select_all_results", "Select All", class = "btn-info", width = "100%")),
                column(4, actionButton("deselect_all_results", "Deselect All", class = "btn-default", width = "100%")),
                column(4, actionButton("select_for_analysis", "Send to Analysis", class = "btn-success", width = "100%"))
              ),
              p(style = "margin-top: 10px;", textOutput("selection_count"))
          )
        )
      ),
      
      # Add Data Tab
      tabItem(tabName = "add",
        fluidRow(
          box(title = "Add Single Sequence", status = "success", solidHeader = TRUE, width = 6,
              textInput("new_sample_id", "Sample ID:", placeholder = "e.g., SAMPLE_001"),
              textAreaInput("new_sequence", "Sequence:", rows = 5, 
                           placeholder = "Paste DNA sequence here..."),
              hr(),
              h4("Metadata (optional)"),
              dateInput("new_date", "Sampling Date:", value = NA),
              textInput("new_genotype", "Genotype:", placeholder = "e.g., IA, IB, IIIA"),
              textInput("new_variant", "Variant:", placeholder = "e.g., Specific variant"),
              textInput("new_patient_id", "Patient ID:", placeholder = "e.g., Patient identifier"),
              textInput("new_location", "Location:", placeholder = "e.g., Oslo"),
              textInput("new_country", "Country:", placeholder = "e.g., Norway"),
              textInput("new_source", "Source:", placeholder = "e.g., Blood, Stool"),
              textInput("new_transmission", "Transmission Route:", placeholder = "e.g., Travel, Food"),
              textAreaInput("new_comment", "Comment:", rows = 2),
              actionButton("add_sequence_btn", "Add Sequence", class = "btn-success", width = "100%")
          ),
          box(title = "Import from FASTA", status = "info", solidHeader = TRUE, width = 6,
              fileInput("fasta_file", "Select FASTA File:", accept = c(".fasta", ".fa", ".fna")),
              checkboxInput("skip_existing", "Skip existing samples", value = TRUE),
              hr(),
              h4("Optional: Import Metadata CSV"),
              fileInput("metadata_file", "Select Metadata CSV:", accept = ".csv"),
              p("CSV must contain 'sample_id' column. Other columns: sampling_date, sample_year, genotype, 
                 variant, patient_id, geo_location, geo_country, source, transmission_route, comment"),
              hr(),
              actionButton("import_btn", "Import", class = "btn-primary", width = "100%"),
              verbatimTextOutput("import_result")
          )
        )
      ),
      
      # Analysis Tab
      tabItem(tabName = "analysis",
        fluidRow(
          box(title = "Selected Sequences for Analysis", status = "primary", solidHeader = TRUE, width = 12,
              DT::dataTableOutput("analysis_sequences"),
              actionButton("clear_analysis_selection", "Clear Selection", class = "btn-warning")
          )
        ),
        fluidRow(
          box(title = "Multiple Sequence Alignment", status = "info", solidHeader = TRUE, width = 6,
              selectInput("msa_method", "Alignment Method:",
                         choices = c("Muscle", "ClustalW", "ClustalOmega")),
              actionButton("run_msa_btn", "Run Alignment", class = "btn-primary"),
              hr(),
              verbatimTextOutput("msa_stats"),
              downloadButton("download_alignment", "Download Alignment (FASTA)")
          ),
          box(title = "Phylogenetic Tree", status = "success", solidHeader = TRUE, width = 6,
              selectInput("tree_method", "Tree Method:",
                         choices = c("Neighbor-Joining" = "nj", "UPGMA" = "upgma")),
              selectInput("dist_model", "Distance Model:",
                         choices = c("K80", "K81", "F81", "F84", "T92", "TN93", "JC69", "raw")),
              actionButton("run_tree_btn", "Build Tree", class = "btn-success"),
              hr(),
              plotOutput("tree_plot", height = 400),
              downloadButton("download_tree", "Download Tree (Newick)")
          )
        ),
        fluidRow(
          box(title = "Save / Manage Analyses", status = "primary", solidHeader = TRUE, width = 6,
              textInput("analysis_name", "Analysis name:", placeholder = "Short descriptive name"),
              textAreaInput("analysis_description", "Description:", rows = 3, placeholder = "Optional description"),
              actionButton("save_analysis", "Save analysis", class = "btn-success", width = "100%"),
              br(), br(),
              actionButton("load_selected_analysis", "Load selected analysis", class = "btn-info", width = "48%"),
              actionButton("delete_selected_analysis", "Delete selected analysis", class = "btn-danger", width = "48%")
          ),
          box(title = "Stored Analyses", status = "info", solidHeader = TRUE, width = 6,
              DT::dataTableOutput("stored_analyses"),
              p(style = "font-size: 11px;", em("Select a row then click 'Load selected analysis' to load into the viewer."))
          )
        ),
        fluidRow(
          box(title = "Distance Matrix", status = "warning", solidHeader = TRUE, width = 6,
              actionButton("run_dist_btn", "Calculate Distances", class = "btn-warning"),
              plotlyOutput("dist_heatmap", height = 400)
          ),
          box(title = "Clustering", status = "danger", solidHeader = TRUE, width = 6,
              numericInput("cluster_threshold", "Distance Threshold:", value = 0.03, 
                          min = 0, max = 1, step = 0.01),
              actionButton("run_cluster_btn", "Cluster Sequences", class = "btn-danger"),
              hr(),
              DT::dataTableOutput("cluster_results")
          )
        )
      ),
      
      # Export Tab
      tabItem(tabName = "export",
        fluidRow(
          box(title = "Export Sequences", status = "primary", solidHeader = TRUE, width = 6,
              p("Export all sequences or selected sequences from search/analysis."),
              checkboxInput("export_include_meta", "Include metadata in FASTA headers", value = FALSE),
              downloadButton("export_fasta", "Export All as FASTA", class = "btn-primary"),
              downloadButton("export_selected_fasta", "Export Selected as FASTA", class = "btn-info")
          ),
          box(title = "Export Metadata", status = "info", solidHeader = TRUE, width = 6,
              p("Export metadata as CSV file."),
              downloadButton("export_csv", "Export All Metadata (CSV)", class = "btn-primary"),
              downloadButton("export_selected_csv", "Export Selected Metadata (CSV)", class = "btn-info")
          )
        )
      ),
      
      # Microreact Tab
      tabItem(tabName = "microreact",
        fluidRow(
          box(title = "Send to Microreact", status = "primary", solidHeader = TRUE, width = 6,
              p("Create a Microreact visualization from your phylogenetic tree and metadata."),
              p(strong("Requirements:"), "Run alignment and build a tree first in the Analysis tab."),
              hr(),
              textInput("microreact_project_name", "Project Name:", 
                       value = paste0("HAV_Analysis_", format(Sys.Date(), "%Y%m%d"))),
              textAreaInput("microreact_description", "Description:", 
                           rows = 2, value = "HAV phylogenetic analysis"),
              hr(),
              actionButton("send_to_microreact", "Send to Microreact", 
                          class = "btn-success", icon = icon("upload"), width = "100%"),
              br(), br(),
              h5("Manual Upload Options:"),
              downloadButton("download_microreact_file", "Download Files (ZIP)", 
                            class = "btn-info", style = "width: 100%;"),
              p(em("Downloads metadata.csv + tree.nwk - upload both to microreact.org/upload"), style = "font-size: 11px; margin-top: 5px;"),
              fluidRow(
                column(6, downloadButton("download_microreact_csv", "Download CSV", class = "btn-default", style = "width: 100%;")),
                column(6, downloadButton("download_microreact_tree", "Download Tree", class = "btn-default", style = "width: 100%;"))
              ),
              hr(),
                actionButton("save_open_local", "Save & Open in Local Viewer", class = "btn-warning", width = "100%"),
                p(em("Writes .microreact into viewer/public/data and opens embedded viewer"), style = "font-size:11px; margin-top:5px;"),
              verbatimTextOutput("microreact_status")
          ),
          box(title = "Microreact Visualization", status = "success", solidHeader = TRUE, width = 6,
              p("After sending to Microreact, the visualization will appear below."),
              uiOutput("microreact_link"),
              hr(),
              uiOutput("microreact_iframe")
          )
        ),
        fluidRow(
          box(title = "Previous Microreact Projects", status = "info", solidHeader = TRUE, width = 12,
              p("Your recent Microreact projects from this session:"),
              DT::dataTableOutput("microreact_history")
          )
        )
      ),
      
      # Settings Tab
      tabItem(tabName = "settings",
        fluidRow(
          box(title = "Microreact API Settings", status = "warning", solidHeader = TRUE, width = 6,
              p("Enter your Microreact API token to enable direct upload."),
              p("Get your token from: ", 
                tags$a(href = "https://microreact.org/my-account/settings", 
                       target = "_blank", "https://microreact.org/my-account/settings")),
              hr(),
              passwordInput("microreact_api_token", "API Token:", value = ""),
              actionButton("save_api_token", "Save Token", class = "btn-primary"),
              actionButton("test_api_token", "Test Connection", class = "btn-info"),
              hr(),
              verbatimTextOutput("api_token_status")
          ),
          box(title = "About Microreact", status = "info", solidHeader = TRUE, width = 6,
              p("Microreact is a web application for visualization of phylogenetic trees with metadata."),
              p("Features:"),
              tags$ul(
                tags$li("Interactive phylogenetic tree visualization"),
                tags$li("Geographic mapping of samples"),
                tags$li("Timeline visualization"),
                tags$li("Metadata table with filtering")
              ),
              p("Learn more: ", 
                tags$a(href = "https://microreact.org", target = "_blank", "https://microreact.org"))
          )
        )
      )
    )
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  
  # Reactive values
  rv <- reactiveValues(
    selected_for_analysis = character(),
    msa_result = NULL,
    tree_result = NULL,
    dist_matrix = NULL,
    microreact_token = NULL,
    microreact_projects = data.frame(
      name = character(),
      url = character(),
      created = character(),
      stringsAsFactors = FALSE
    ),
    current_microreact_url = NULL
    ,
    auto_open_local = FALSE
  )
  
  # Database connection (per session)
  con <- reactive({
    if (!file.exists(DB_PATH)) {
      showNotification("Database not found. Please run setup_database.R first.", type = "error")
      return(NULL)
    }
    get_db_connection(DB_PATH)
  })
  
  # Clean up on session end
  session$onSessionEnded(function() {
    if (!is.null(con())) {
      try(close_db_connection(con()), silent = TRUE)
    }
  })
  
  # ============================================
  # Dashboard
  # ============================================
  
  output$total_sequences <- renderInfoBox({
    req(con())
    n <- dbGetQuery(con(), "SELECT COUNT(*) as n FROM sequences")$n
    infoBox("Total Sequences", n, icon = icon("dna"), color = "blue")
  })
  
  output$total_viruses <- renderInfoBox({
    req(con())
    n <- dbGetQuery(con(), "SELECT COUNT(DISTINCT genotype) as n FROM metadata WHERE genotype IS NOT NULL")$n
    infoBox("Genotypes", n, icon = icon("virus"), color = "green")
  })

  output$date_range <- renderInfoBox({
    req(con())
    years <- dbGetQuery(con(), "SELECT MIN(sample_year) as min_y, MAX(sample_year) as max_y FROM metadata")
    range_str <- if (!is.na(years$min_y)) paste(years$min_y, "-", years$max_y) else "No years"
    infoBox("Year Range", range_str, icon = icon("calendar"), color = "purple")
  })
  
  output$virus_plot <- renderPlotly({
    req(con())
    data <- dbGetQuery(con(), "SELECT genotype, COUNT(*) as n FROM metadata WHERE genotype IS NOT NULL GROUP BY genotype")
    if (nrow(data) == 0) return(NULL)
    plot_ly(data, x = ~genotype, y = ~n, type = "bar") %>%
      layout(xaxis = list(title = "Genotype"), yaxis = list(title = "Count"))
  })
  
  output$location_plot <- renderPlotly({
    req(con())
    data <- dbGetQuery(con(), "SELECT geo_country, COUNT(*) as n FROM metadata WHERE geo_country IS NOT NULL GROUP BY geo_country ORDER BY n DESC LIMIT 20")
    if (nrow(data) == 0) return(NULL)
    plot_ly(data, x = ~geo_country, y = ~n, type = "bar", marker = list(color = "green")) %>%
      layout(xaxis = list(title = "Country", tickangle = 45), yaxis = list(title = "Count"))
  })
  
  output$recent_sequences <- renderDT({
    req(con())
    data <- dbGetQuery(con(), "
      SELECT s.sample_id, s.sequence_length, m.genotype, m.variant, m.geo_country, m.sample_year
      FROM sequences s LEFT JOIN metadata m ON s.sample_id = m.sample_id
      ORDER BY s.created_at DESC LIMIT 10
    ")
    datatable(data, options = list(pageLength = 5, searching = FALSE), rownames = FALSE)
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
  
  # ============================================
  # Analysis
  # ============================================
  
  output$analysis_sequences <- renderDT({
    if (length(rv$selected_for_analysis) == 0) {
      return(data.frame(Message = "No sequences selected. Use Search tab to select sequences."))
    }
    req(con())
    data <- get_sequences_with_metadata(con(), rv$selected_for_analysis)
    data$sequence <- NULL  # Don't show sequence
    datatable(data, selection = "none", options = list(pageLength = 10), rownames = FALSE)
  })
  
  observeEvent(input$clear_analysis_selection, {
    rv$selected_for_analysis <- character()
    rv$msa_result <- NULL
    rv$tree_result <- NULL
    rv$dist_matrix <- NULL
  })
  
  # MSA
  observeEvent(input$run_msa_btn, {
    req(length(rv$selected_for_analysis) >= 2)
    req(con())
    
    withProgress(message = "Running alignment...", {
      tryCatch({
        rv$msa_result <- run_msa_from_db(con(), rv$selected_for_analysis, input$msa_method)
        showNotification("Alignment complete!", type = "message")

      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
  
  output$msa_stats <- renderPrint({
    req(rv$msa_result)
    stats <- alignment_stats(rv$msa_result)
    cat("Alignment Statistics\n")
    cat("====================\n")
    cat("Sequences:", stats$n_sequences, "\n")
    cat("Alignment length:", stats$alignment_length, "bp\n")
    cat("Conserved positions:", stats$conserved_positions, "(", stats$conservation_pct, "%)\n")
    cat("Gap positions:", stats$gap_positions, "\n")
  })
  
  output$download_alignment <- downloadHandler(
    filename = function() paste0("alignment_", Sys.Date(), ".fasta"),
    content = function(file) {
      req(rv$msa_result)
      export_alignment_fasta(rv$msa_result, file)
    }
  )
  
  # Tree
  observeEvent(input$run_tree_btn, {
    req(rv$msa_result)
    
    withProgress(message = "Building tree...", {
      tryCatch({
        if (input$tree_method == "nj") {
          rv$tree_result <- build_nj_tree(rv$msa_result, input$dist_model)
        } else {
          rv$tree_result <- build_upgma_tree(rv$msa_result, input$dist_model)
        }
        showNotification("Tree built!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
  
  output$tree_plot <- renderPlot({
    req(rv$tree_result)
    plot(rv$tree_result, type = "phylogram", cex = 0.8)
    title(paste(input$tree_method, "tree"))
  })
  
  output$download_tree <- downloadHandler(
    filename = function() paste0("tree_", Sys.Date(), ".nwk"),
    content = function(file) {
      req(rv$tree_result)
      export_tree_newick(rv$tree_result, file)
    }
  )
  
  # Distance matrix
  observeEvent(input$run_dist_btn, {
    req(rv$msa_result)
    
    tryCatch({
      rv$dist_matrix <- calculate_distance_matrix(rv$msa_result, input$dist_model)
      showNotification("Distance matrix calculated!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$dist_heatmap <- renderPlotly({
    req(rv$dist_matrix)
    mat <- as.matrix(rv$dist_matrix)
    plot_ly(z = mat, x = colnames(mat), y = rownames(mat), type = "heatmap",
            colors = colorRamp(c("white", "blue"))) %>%
      layout(xaxis = list(tickangle = 45))
  })
  
  # Clustering
  observeEvent(input$run_cluster_btn, {
    req(rv$msa_result)
    
    tryCatch({
      clusters <- cluster_by_distance(rv$msa_result, input$cluster_threshold, input$dist_model)
      
      output$cluster_results <- renderDT({
        datatable(clusters, options = list(pageLength = 10), rownames = FALSE)
      })
      
      showNotification(paste("Found", length(unique(clusters$cluster)), "clusters"), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
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
  
  # ============================================
  # Microreact Settings
  # ============================================
  
  # Save API token
  observeEvent(input$save_api_token, {
    if (nchar(input$microreact_api_token) > 0) {
      rv$microreact_token <- input$microreact_api_token
      showNotification("API token saved for this session", type = "message")
      output$api_token_status <- renderPrint({
        cat("✓ Token saved\n")
        cat("Token will be validated when creating a project.\n")
      })
    } else {
      showNotification("Please enter an API token", type = "warning")
    }
  })
  
  # Test API token - validate format and save
  observeEvent(input$test_api_token, {
    req(input$microreact_api_token)
    
    token <- input$microreact_api_token
    
    # Check if token looks like a JWT (starts with eyJ)
    if (grepl("^eyJ", token) && nchar(token) > 50) {
      rv$microreact_token <- token
      output$api_token_status <- renderPrint({
        cat("✓ Token format looks valid (JWT)\n")
        cat("Token saved for this session.\n")
        cat("Note: Full validation occurs when creating a project.\n")
        cat("\nToken preview:", substr(token, 1, 20), "...\n")
      })
      showNotification("Token saved! Will be validated on first use.", type = "message")
    } else {
      output$api_token_status <- renderPrint({
        cat("⚠ Token format may be incorrect\n")
        cat("Microreact tokens typically start with 'eyJ' (JWT format)\n")
        cat("Please verify you copied the complete token from:\n")
        cat("https://microreact.org/my-account/settings\n")
      })
      showNotification("Token format looks unusual - please verify", type = "warning")
    }
  })
  
  # ============================================
  # Microreact Integration
  # ============================================
  
  # Helper function to create microreact project data
  create_microreact_data <- function() {
    req(rv$tree_result)
    req(con())
    
    # Get metadata for selected samples
    sample_ids <- rv$selected_for_analysis
    metadata <- get_sequences_with_metadata(con(), sample_ids)
    
    # Get tree tip labels - these are what will be in the Newick file
    # write.tree() converts spaces to underscores, so we need to match
    tree_tip_labels <- rv$tree_result$tip.label
    
    # Create a mapping from original sample_id to tree tip label
    # The tree was built from sequences, so tip labels match the MSA row names
    # which are the sample_ids with spaces converted to underscores
    id_to_tree_label <- function(id) {
      gsub(" ", "_", id)
    }
    
    # Prepare metadata CSV for Microreact
    # Use the tree tip label format for the id column so it matches the Newick file
    microreact_meta <- data.frame(
      id = sapply(metadata$sample_id, id_to_tree_label),
      original_id = metadata$sample_id,
      genotype = ifelse(is.na(metadata$genotype), "", metadata$genotype),
      variant = ifelse(is.na(metadata$variant), "", metadata$variant),
      patient_id = ifelse(is.na(metadata$patient_id), "", metadata$patient_id),
      country = ifelse(is.na(metadata$geo_country), "", metadata$geo_country),
      location = ifelse(is.na(metadata$geo_location), "", metadata$geo_location),
      year = ifelse(is.na(metadata$sample_year), "", as.character(metadata$sample_year)),
      stringsAsFactors = FALSE
    )

    # Deterministic country colour mapping (adds country__colour column)
    countries <- unique(na.omit(microreact_meta$country[microreact_meta$country != ""]))
    if (length(countries) > 0) {
      n <- length(countries)
      if (requireNamespace("RColorBrewer", quietly = TRUE)) {
        pal_size <- min(12, max(3, n))
        pal <- RColorBrewer::brewer.pal(pal_size, "Set3")
        cols <- rep_len(pal, n)
      } else {
        cols <- grDevices::rainbow(n)
      }
      country_col_map <- setNames(cols, countries)
      microreact_meta$country__colour <- vapply(microreact_meta$country, function(x) {
        if (is.na(x) || x == "") return("")
        if (!is.null(country_col_map[[x]])) return(country_col_map[[x]])
        return("")
      }, FUN.VALUE = "")
    } else {
      microreact_meta$country__colour <- ""
    }

    # Try to normalize latitude / longitude columns if present in source metadata
    find_col <- function(patterns, df) {
      nm <- names(df)
      for (p in patterns) {
        idx <- which(tolower(nm) == tolower(p))
        if (length(idx)) return(nm[idx[1]])
      }
      for (p in patterns) {
        idx <- grep(p, tolower(nm))
        if (length(idx)) return(nm[idx[1]])
      }
      return(NA)
    }

    lat_col <- find_col(c("latitude", "lat", "y", "geo_lat", "gps_lat"), metadata)
    lon_col <- find_col(c("longitude", "lon", "lng", "x", "geo_lon", "gps_lon"), metadata)

    if (!is.na(lat_col) && !is.na(lon_col)) {
      microreact_meta$latitude <- as.numeric(metadata[[lat_col]])
      microreact_meta$longitude <- as.numeric(metadata[[lon_col]])
    } else {
      microreact_meta$latitude <- NA
      microreact_meta$longitude <- NA
    }
    
    # Convert tree to Newick string
    tree_newick <- write.tree(rv$tree_result)
    
    # Convert metadata to CSV string
    temp_csv <- tempfile(fileext = ".csv")
    write.csv(microreact_meta, temp_csv, row.names = FALSE, na = "")
    meta_csv <- paste(readLines(temp_csv, warn = FALSE), collapse = "\n")
    unlink(temp_csv)
    
    # Return as a list with separate components for flexibility
    list(
      metadata_csv = meta_csv,
      tree_newick = tree_newick,
      metadata_df = microreact_meta
    )
  }
  
  # Download as separate files (ZIP with CSV + Newick)
  output$download_microreact_file <- downloadHandler(
    filename = function() {
      paste0(input$microreact_project_name, "_microreact.zip")
    },
    content = function(file) {
      tryCatch({
        data <- create_microreact_data()
        
        # Create temp directory for files
        temp_dir <- tempdir()
        csv_file <- file.path(temp_dir, "metadata.csv")
        tree_file <- file.path(temp_dir, "tree.nwk")
        
        # Write files
        writeLines(data$metadata_csv, csv_file)
        writeLines(data$tree_newick, tree_file)
        
        # Create ZIP
        old_wd <- getwd()
        setwd(temp_dir)
        zip(file, files = c("metadata.csv", "tree.nwk"))
        setwd(old_wd)
        
        showNotification("Files downloaded! Extract ZIP and upload both files to microreact.org/upload", type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    }
  )
  
  # Download CSV file separately for Microreact
  output$download_microreact_csv <- downloadHandler(
    filename = function() {
      paste0(input$microreact_project_name, "_metadata.csv")
    },
    content = function(file) {
      tryCatch({
        data <- create_microreact_data()
        writeLines(data$metadata_csv, file)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    }
  )
  
  # Download Newick tree file separately for Microreact
  output$download_microreact_tree <- downloadHandler(
    filename = function() {
      paste0(input$microreact_project_name, "_tree.nwk")
    },
    content = function(file) {
      tryCatch({
        data <- create_microreact_data()
        writeLines(data$tree_newick, file)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    }
  )
  
  # Send to Microreact
  observeEvent(input$send_to_microreact, {
    req(rv$tree_result)
    req(rv$microreact_token)
    req(con())
    
    tryCatch({
      showNotification("Preparing data for Microreact...", type = "message", duration = 2)
      
      # Use the helper function to get properly formatted data
      data <- create_microreact_data()
      
      # Generate unique IDs for files
      timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
      data_file_id <- paste0("data-", timestamp)
      tree_file_id <- paste0("tree-", timestamp)
      
      # Create proper .microreact JSON structure for API
      files_obj <- list()
      files_obj[[data_file_id]] <- list(name = "metadata.csv", format = "text/csv", blob = data$metadata_csv)
      files_obj[[tree_file_id]] <- list(name = "tree.nwk", format = "text/x-nh", blob = data$tree_newick)

      microreact_project <- list(
        meta = list(
          name = input$microreact_project_name,
          description = input$microreact_description
        ),
        files = files_obj,
        datasets = list(
          list(
            id = "dataset-1",
            file = data_file_id,
            idFieldName = "id"
          )
        ),
        trees = list(
          list(
            id = "tree-1",
            file = tree_file_id,
            labelField = "id"
          )
        )
        ,
        # Provide viewer-friendly initial view preferences: label and colour by country
        initial_view = list(
          labelField = "country",
          colourField = "country",
          colourColumn = "country__colour",
          showMap = TRUE
        )
      )
      
      # Convert to JSON
      json_body <- toJSON(microreact_project, auto_unbox = TRUE)
      
      showNotification("Uploading to Microreact...", type = "message", duration = 2)
      
      # Send to Microreact API
      response <- POST(
        url = "https://microreact.org/api/projects/create",
        add_headers(
          `Access-Token` = rv$microreact_token,
          `Content-Type` = "application/json; charset=utf-8"
        ),
        body = json_body,
        encode = "raw"
      )
      
      if (status_code(response) == 200) {
        result <- content(response, as = "parsed")
        project_url <- result$url
        project_id <- result$id

        # Store the URL
        rv$current_microreact_url <- project_url

        # Add to history
        rv$microreact_projects <- rbind(
          data.frame(
            name = input$microreact_project_name,
            url = project_url,
            created = format(Sys.time(), "%Y-%m-%d %H:%M"),
            stringsAsFactors = FALSE
          ),
          rv$microreact_projects
        )
        showNotification("Microreact project created!", type = "message")
      } else {
        # Get detailed error info
        error_content <- content(response, as = "text")
        
        output$microreact_status <- renderPrint({
          cat("✗ Failed to create project\n")
          cat("Status:", status_code(response), "\n")
          cat("Response:", error_content, "\n\n")
          cat("TIP: Try downloading the .microreact file and uploading manually.\n")
        })
        showNotification("Failed to create Microreact project", type = "error")
      }
      
    }, error = function(e) {
      output$microreact_status <- renderPrint({
        cat("✗ Error\n")
        cat(e$message, "\n")
      })
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Render Microreact link
  output$microreact_link <- renderUI({
    if (!is.null(rv$current_microreact_url)) {
      tagList(
        p(strong("Project URL:")),
        tags$a(href = rv$current_microreact_url, target = "_blank", rv$current_microreact_url),
        br(),
        actionButton("open_microreact_iframe", "Load in Viewer Below", class = "btn-info btn-sm")
      )
    } else {
      p("No project created yet. Run analysis and click 'Send to Microreact'.")
    }
  })
  
  # Render Microreact iframe
  output$microreact_iframe <- renderUI({
    if (!is.null(rv$current_microreact_url) && ((!is.null(input$open_microreact_iframe) && input$open_microreact_iframe > 0) || isTRUE(rv$auto_open_local))) {
      tags$iframe(
        src = rv$current_microreact_url,
        width = "100%",
        height = "600px",
        frameborder = "0",
        style = "border: 1px solid #ddd; border-radius: 4px;"
      )
    } else if (!is.null(rv$current_microreact_url)) {
      p("Click 'Load in Viewer Below' to display the visualization, or open the link in a new tab.")
    } else {
      p(em("Visualization will appear here after creating a project."))
    }
  })

  # Save .microreact into local viewer public/data and open in embedded viewer
  observeEvent(input$save_open_local, {
    req(rv$tree_result)
    tryCatch({
      data <- create_microreact_data()

      # Build .microreact JSON
      timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
      safe_name <- gsub("[^A-Za-z0-9._-]", "_", input$microreact_project_name)
      filename <- paste0(safe_name, "_", timestamp, ".microreact")

      data_file_id <- paste0("data-", timestamp)
      tree_file_id <- paste0("tree-", timestamp)

      files_obj <- list()
      files_obj[[data_file_id]] <- list(name = "metadata.csv", format = "text/csv", blob = data$metadata_csv)
      files_obj[[tree_file_id]] <- list(name = "tree.nwk", format = "text/x-nh", blob = data$tree_newick)

      microreact_project <- list(
        meta = list(
          name = input$microreact_project_name,
          description = input$microreact_description
        ),
        files = files_obj,
        datasets = list(list(id = "dataset-1", file = data_file_id, idFieldName = "id")),
        trees = list(list(id = "tree-1", file = tree_file_id, labelField = "id"))
        ,
        initial_view = list(
          labelField = "country",
          colourField = "country",
          colourColumn = "country__colour",
          showMap = TRUE
        )
      )

      json_body <- toJSON(microreact_project, auto_unbox = TRUE)

      # Ensure viewer public data dir exists
      viewer_data_dir <- file.path(app_dir, "viewer", "public", "data")
      dir.create(viewer_data_dir, recursive = TRUE, showWarnings = FALSE)
      out_path <- file.path(viewer_data_dir, filename)
      writeLines(json_body, out_path)

      # Set URL to local viewer with file param
      rv$current_microreact_url <- paste0("http://localhost:3000/?file=", URLencode(filename))
      rv$auto_open_local <- TRUE

      showNotification(paste("Saved .microreact to", out_path, "and set viewer URL."), type = "message")
    }, error = function(e) {
      showNotification(paste("Error saving local .microreact:", e$message), type = "error")
    })
  })
  
  # Microreact history table
  output$microreact_history <- renderDT({
    if (nrow(rv$microreact_projects) > 0) {
      df <- rv$microreact_projects
      df$url <- paste0('<a href="', df$url, '" target="_blank">Open</a>')
      datatable(df, escape = FALSE, options = list(pageLength = 5), rownames = FALSE)
    } else {
      datatable(data.frame(Message = "No projects created in this session"), 
                options = list(dom = 't'), rownames = FALSE)
    }
  })

  # --------------------------------------------
  # Analysis persistence (save / load analyses)
  # --------------------------------------------

  analysis_dir <- file.path(app_dir, "analysis")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

  list_saved_analyses <- function() {
    files <- list.files(analysis_dir, pattern = "\\.rds$", full.names = TRUE)
    infos <- lapply(files, function(f) {
      ok <- tryCatch({ a <- readRDS(f); TRUE }, error = function(e) FALSE)
      if (!ok) return(NULL)
      a <- readRDS(f)
      data.frame(
        file = basename(f),
        name = if (!is.null(a$name)) a$name else "",
        type = if (!is.null(a$type)) a$type else "",
        created = if (!is.null(a$timestamp)) as.character(a$timestamp) else file.info(f)$ctime,
        description = if (!is.null(a$description)) a$description else "",
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, Filter(Negate(is.null), infos))
  }

  output$stored_analyses <- renderDT({
    df <- list_saved_analyses()
    if (is.null(df) || nrow(df) == 0) return(datatable(data.frame(Message = "No saved analyses"), options = list(dom = 't')))
    datatable(df, selection = 'single', rownames = FALSE, options = list(pageLength = 10))
  })

  observeEvent(input$save_analysis, {
    # Save current analysis (msa and/or tree)
    req(length(rv$selected_for_analysis) > 0)
    name <- ifelse(nzchar(input$analysis_name), input$analysis_name, paste0("analysis_", format(Sys.time(), "%Y%m%d%H%M%S")))
    description <- input$analysis_description %||% ""
    a_type <- if (!is.null(rv$tree_result)) "phylogeny" else if (!is.null(rv$msa_result)) "alignment" else "analysis"
    obj <- list(
      name = name,
      description = description,
      type = a_type,
      timestamp = Sys.time(),
      sample_ids = rv$selected_for_analysis,
      msa = rv$msa_result,
      tree = rv$tree_result
    )
    safe <- gsub("[^A-Za-z0-9._-]", "_", name)
    filename <- paste0(safe, "_", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")
    outpath <- file.path(analysis_dir, filename)
    saveRDS(obj, outpath)
    showNotification(paste("Saved analysis:", filename), type = "message")
    # trigger refresh
    output$stored_analyses <- renderDT({
      df <- list_saved_analyses()
      datatable(df, selection = 'single', rownames = FALSE, options = list(pageLength = 10))
    })
  })

  observeEvent(input$load_selected_analysis, {
    sel <- input$stored_analyses_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showNotification("No analysis selected", type = "warning")
      return()
    }
    df <- list_saved_analyses()
    row <- df[sel, ]
    if (is.null(row) || nrow(row) == 0) {
      showNotification("Selection missing", type = "error")
      return()
    }
    path <- file.path(analysis_dir, row$file)
    obj <- tryCatch(readRDS(path), error = function(e) { NULL })
    if (is.null(obj)) {
      showNotification("Failed to read analysis file", type = "error")
      return()
    }
    # Load into reactive values
    rv$msa_result <- obj$msa
    rv$tree_result <- obj$tree
    rv$selected_for_analysis <- obj$sample_ids %||% rv$selected_for_analysis
    showNotification(paste("Loaded analysis:", obj$name), type = "message")
  })

  observeEvent(input$delete_selected_analysis, {
    sel <- input$stored_analyses_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showNotification("No analysis selected", type = "warning")
      return()
    }
    df <- list_saved_analyses()
    row <- df[sel, ]
    path <- file.path(analysis_dir, row$file)
    if (file.exists(path)) {
      file.remove(path)
      showNotification(paste("Deleted:", row$file), type = "message")
    }
    # refresh
    output$stored_analyses <- renderDT({
      df <- list_saved_analyses()
      if (is.null(df) || nrow(df) == 0) return(datatable(data.frame(Message = "No saved analyses"), options = list(dom = 't')))
      datatable(df, selection = 'single', rownames = FALSE, options = list(pageLength = 10))
    })
  })
}

# =============================================================================
# Run App
# =============================================================================

shinyApp(ui = ui, server = server)
