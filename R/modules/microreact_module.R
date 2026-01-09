register_microreact_server <- function(input, output, session, rv, con, app_dir){
  # Source the API functions
  source(file.path(app_dir, "R", "microreact_api.R"), local = TRUE)
  
  # File path for persistent token storage (not in git)
  token_file <- file.path(app_dir, ".microreact_token")
  
  # Load saved token on startup
  if (file.exists(token_file)) {
    saved_token <- tryCatch({
      readLines(token_file, warn = FALSE)[1]
    }, error = function(e) NULL)
    if (!is.null(saved_token) && nchar(saved_token) > 50) {
      rv$microreact_token <- saved_token
      # Update the input field with saved token
      updateTextInput(session, "microreact_api_token", value = saved_token)
      output$api_token_status <- renderPrint({
        cat("✓ Token loaded from saved settings\n")
        cat("Token preview:", substr(saved_token, 1, 20), "...\n")
      })
    }
  }
  
  # Microreact Settings handlers
  observeEvent(input$save_api_token, {
    if (nchar(input$microreact_api_token) > 0) {
      rv$microreact_token <- input$microreact_api_token
      # Persist to file
      tryCatch({
        writeLines(input$microreact_api_token, token_file)
        showNotification("API token saved permanently", type = "message")
        output$api_token_status <- renderPrint({
          cat("✓ Token saved permanently\n")
          cat("Token will be loaded automatically on restart.\n")
        })
      }, error = function(e) {
        showNotification("Token saved for session only (could not persist to file)", type = "warning")
        output$api_token_status <- renderPrint({
          cat("✓ Token saved for this session\n")
          cat("⚠ Could not persist to file:", e$message, "\n")
        })
      })
    } else {
      showNotification("Please enter an API token", type = "warning")
    }
  })

  observeEvent(input$test_api_token, {
    req(input$microreact_api_token)
    token <- input$microreact_api_token
    if (grepl("^eyJ", token) && nchar(token) > 50) {
      rv$microreact_token <- token
      # Also persist when testing
      tryCatch({
        writeLines(token, token_file)
      }, error = function(e) NULL)
      output$api_token_status <- renderPrint({
        cat("✓ Token format looks valid (JWT)\n")
        cat("Token saved permanently.\n")
        cat("Note: Full validation occurs when creating a project.\n")
        cat("\nToken preview:", substr(token, 1, 20), "...\n")
      })
      showNotification("Token saved permanently!", type = "message")
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

  # Helper to build microreact data
  create_microreact_data <- function(){
    req(rv$tree_result)
    req(con())
    sample_ids <- rv$selected_for_analysis
    metadata <- get_sequences_with_metadata(con(), sample_ids)
    tree_tip_labels <- rv$tree_result$tip.label
    id_to_tree_label <- function(id) { gsub(" ", "_", id) }
    
    # Load country codes from JSON file (contains country name, alpha2, alpha3, lat, long)
    country_codes_file <- file.path(app_dir, "country-codes-lat-long-alpha3.json")
    country_lookup <- NULL
    if (file.exists(country_codes_file)) {
      tryCatch({
        json_data <- fromJSON(country_codes_file)
        country_lookup <- json_data$ref_country_codes
      }, error = function(e) {
        warning("Could not load country codes JSON: ", e$message)
      })
    }
    
    # Additional country name aliases (maps local names to standard names in JSON)
    country_aliases <- c(
      "Marokko" = "Morocco",
      "Sverige" = "Sweden", 
      "Tyskland" = "Germany",
      "Türkiye" = "Turkey",
      "Czech Republic" = "Czechia",
      "UK" = "United Kingdom",
      "USA" = "United States",
      "UAE" = "United Arab Emirates",
      "South Korea" = "Korea, Republic of",
      "Korea" = "Korea, Republic of",
      "Russia" = "Russian Federation",
      "Iran" = "Iran, Islamic Republic of",
      "Syria" = "Syrian Arab Republic",
      "Vietnam" = "Viet Nam",
      "Tanzania" = "Tanzania, United Republic of",
      "Scotland" = "United Kingdom",
      "England" = "United Kingdom",
      "Wales" = "United Kingdom"
    )
    
    # Function to lookup country info (returns list with alpha2, latitude, longitude)
    get_country_info <- function(country_name) {
      result <- list(alpha2 = "", latitude = NA, longitude = NA)
      if (is.na(country_name) || country_name == "" || country_name == "Unknown") return(result)
      
      # Check if this is an alias first
      lookup_name <- country_name
      if (country_name %in% names(country_aliases)) {
        lookup_name <- country_aliases[[country_name]]
      }
      
      if (!is.null(country_lookup)) {
        # Try exact match first
        idx <- which(tolower(country_lookup$country) == tolower(lookup_name))
        if (length(idx) == 0) {
          # Try partial match
          idx <- which(sapply(country_lookup$country, function(x) grepl(lookup_name, x, ignore.case = TRUE)))
        }
        if (length(idx) == 0) {
          # Try reverse partial match (lookup name contains country from JSON)
          idx <- which(sapply(country_lookup$country, function(x) grepl(x, lookup_name, ignore.case = TRUE)))
        }
        
        if (length(idx) > 0) {
          match <- country_lookup[idx[1], ]
          result$alpha2 <- match$alpha2
          result$latitude <- match$latitude
          result$longitude <- match$longitude
        }
      }
      
      return(result)
    }
    
    microreact_meta <- data.frame(
      id = sapply(metadata$sample_id, id_to_tree_label),
      genotype = ifelse(is.na(metadata$genotype), "", metadata$genotype),
      variant = ifelse(is.na(metadata$variant), "", metadata$variant),
      patient_id = ifelse(is.na(metadata$patient_id), "", metadata$patient_id),
      country = ifelse(is.na(metadata$geo_country), "", metadata$geo_country),
      location = ifelse(is.na(metadata$geo_location), "", metadata$geo_location),
      year = ifelse(is.na(metadata$sample_year), "", as.character(metadata$sample_year)),
      stringsAsFactors = FALSE
    )
    
    # Lookup country info (ISO code and coordinates) for each sample
    country_info <- lapply(microreact_meta$country, get_country_info)
    
    # Add ISO 3166-1 alpha-2 codes
    iso_codes <- sapply(country_info, function(x) x$alpha2)
    if (any(iso_codes != "")) {
      microreact_meta$iso_country <- iso_codes
    }
    
    # Add latitude and longitude from country lookup
    latitudes <- sapply(country_info, function(x) x$latitude)
    longitudes <- sapply(country_info, function(x) x$longitude)
    if (any(!is.na(latitudes))) {
      microreact_meta$latitude <- latitudes
      microreact_meta$longitude <- longitudes
    }
    
    tree_newick <- write.tree(rv$tree_result)
    temp_csv <- tempfile(fileext = ".csv")
    write.csv(microreact_meta, temp_csv, row.names = FALSE, na = "")
    meta_csv <- paste(readLines(temp_csv, warn = FALSE), collapse = "\n")
    unlink(temp_csv)
    list(metadata_csv = meta_csv, tree_newick = tree_newick, metadata_df = microreact_meta)
  }
  
  # Helper function to build microreact project JSON with map configuration
  build_microreact_project <- function(data, project_name, description) {
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    data_file_id <- paste0("data", timestamp)
    tree_file_id <- paste0("tree", timestamp)
    
    # Base64 encode the file contents with proper data URL prefix
    csv_base64 <- base64enc::base64encode(charToRaw(data$metadata_csv))
    tree_base64 <- base64enc::base64encode(charToRaw(data$tree_newick))
    
    # Build files as an object with file IDs as keys (matching Microreact schema)
    files_obj <- list()
    files_obj[[data_file_id]] <- list(
      id = data_file_id,
      name = "metadata.csv",
      format = "text/csv",
      type = "data",
      size = nchar(data$metadata_csv),
      blob = paste0("data:text/csv;base64,", csv_base64)
    )
    files_obj[[tree_file_id]] <- list(
      id = tree_file_id,
      name = "tree.nwk",
      format = "text/x-nh",
      type = "tree",
      size = nchar(data$tree_newick),
      blob = paste0("data:application/octet-stream;base64,", tree_base64)
    )
    
    # Check what geographic data is available in metadata
    has_iso <- "iso_country" %in% names(data$metadata_df) && any(data$metadata_df$iso_country != "")
    has_latlong <- "latitude" %in% names(data$metadata_df) && "longitude" %in% names(data$metadata_df)
    
    # Build map configuration based on available data
    map_config <- NULL
    if (has_iso || has_latlong) {
      map_config <- list(
        `map-1` = list(
          id = "map-1",
          title = "Map",
          controls = TRUE,
          showMarkers = TRUE,
          showRegions = TRUE,
          showRegionOutlines = TRUE,
          grouped = TRUE,
          nodeSize = 14,
          minNodeSize = 4,
          maxNodeSize = 64,
          markersOpacity = 100,
          scaleMarkers = FALSE,
          scaleType = "logarithmic",
          type = "mapbox",
          style = "",
          coordinateUnit = "decimal-degrees",
          viewport = list(
            longitude = 0,
            latitude = 0,
            zoom = 1.5,
            pitch = 0,
            bearing = 0,
            padding = list(top = 0, bottom = 0, left = 0, right = 0)
          )
        )
      )
      
      # Configure for lat/long (preferred when available) or ISO 3166 codes
      if (has_latlong) {
        map_config$`map-1`$latitudeField <- "latitude"
        map_config$`map-1`$longitudeField <- "longitude"
      }
      if (has_iso) {
        map_config$`map-1`$dataType <- "iso-3166-codes"
        map_config$`map-1`$iso3166Field <- "iso_country"
      }
    }
    
    # Build tree configuration with blocks for country visualization
    tree_config <- list(
      `tree-1` = list(
        id = "tree-1",
        title = "Tree",
        file = tree_file_id,
        labelField = "id",
        type = "rc",
        alignLabels = TRUE,
        showLabels = TRUE,
        showLeafLabels = FALSE,
        showShapes = TRUE,
        showShapeBorders = TRUE,
        showPiecharts = TRUE,
        showEdges = TRUE,
        nodeSize = 14,
        fontSize = 16,
        controls = TRUE,
        blocks = if (has_iso) list("iso_country", "country") else list("country"),
        showBlockHeaders = TRUE,
        blockSize = 14
      )
    )
    
    # Build timeline configuration
    timeline_config <- list(
      `timeline-1` = list(
        id = "timeline-1",
        title = "Timeline",
        dataType = "year-month-day",
        yearField = "year",
        nodeSize = 14,
        style = "bar",
        controls = FALSE
      )
    )
    
    # Build table configuration
    table_config <- list(
      `table-1` = list(
        id = "table-1",
        title = "Metadata",
        file = data_file_id,
        displayMode = "cosy",
        hideUnselected = FALSE
      )
    )
    
    # Build panes layout
    panes_model <- list(
      global = list(
        splitterSize = 2,
        tabEnableClose = FALSE,
        tabSetHeaderHeight = 1,
        tabSetTabStripHeight = 1,
        tabSetMinWidth = 160,
        tabSetMinHeight = 160,
        borderMinSize = 160,
        borderBarSize = 20,
        borderEnableDrop = FALSE
      ),
      borders = list(
        list(
          type = "border",
          size = 240,
          location = "right",
          children = list(
            list(type = "tab", id = "--mr-legend-pane", name = "Legend", component = "Legend", enableClose = FALSE, enableDrag = FALSE),
            list(type = "tab", id = "--mr-selection-pane", name = "Selection", component = "Selection", enableClose = FALSE, enableDrag = FALSE),
            list(type = "tab", id = "--mr-history-pane", name = "History", component = "History", enableClose = FALSE, enableDrag = FALSE),
            list(type = "tab", id = "--mr-views-pane", name = "Views", component = "Views", enableClose = FALSE, enableDrag = FALSE)
          )
        )
      ),
      layout = list(
        type = "row",
        id = "#main-row",
        children = list(
          list(
            type = "row",
            id = "#content-row",
            children = list(
              list(
                type = "row",
                id = "#top-row",
                weight = 64,
                children = list(
                  list(type = "tabset", id = "#map-tabset", children = list(
                    list(type = "tab", id = "map-1", name = "Map", component = "Map")
                  )),
                  list(type = "tabset", id = "#tree-tabset", children = list(
                    list(type = "tab", id = "tree-1", name = "Tree", component = "Tree")
                  ))
                )
              ),
              list(
                type = "tabset",
                id = "#bottom-tabset",
                weight = 48,
                children = list(
                  list(type = "tab", id = "timeline-1", name = "Timeline", component = "Timeline"),
                  list(type = "tab", id = "table-1", name = "Metadata", component = "Table")
                )
              )
            )
          )
        )
      )
    )
    
    microreact_project <- list(
      schema = "https://microreact.org/schema/v1.json",
      meta = list(
        name = project_name, 
        description = description,
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
      ),
      files = files_obj,
      datasets = list(`dataset-1` = list(id = "dataset-1", file = data_file_id, idFieldName = "id")),
      charts = list(),
      filters = list(
        paneFilters = list(),
        dataFilters = list(),
        chartFilters = list(),
        searchOperator = "includes",
        searchValue = "",
        selection = list(),
        selectionBreakdownField = NULL
      ),
      maps = map_config,
      matrices = list(),
      networks = list(),
      notes = list(),
      panes = list(model = panes_model),
      slicers = list(),
      styles = list(
        coloursField = NULL,
        colourPalettes = list(),
        defaultColour = "transparent",
        defaultShape = "circle",
        colourSettings = list(),
        labelsField = NULL,
        legendDirection = "row",
        shapesField = NULL,
        shapePalettes = list()
      ),
      tables = table_config,
      timelines = timeline_config,
      trees = tree_config,
      views = list()
    )
    
    toJSON(microreact_project, auto_unbox = TRUE, null = "null")
  }

  # Downloads: ZIP / CSV / Tree
  output$download_microreact_file <- downloadHandler(
    filename = function() paste0(input$microreact_project_name, "_microreact.zip"),
    content = function(file) {
      tryCatch({
        data <- create_microreact_data()
        temp_dir <- tempdir()
        csv_file <- file.path(temp_dir, "metadata.csv")
        tree_file <- file.path(temp_dir, "tree.nwk")
        writeLines(data$metadata_csv, csv_file)
        writeLines(data$tree_newick, tree_file)
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

  output$download_microreact_csv <- downloadHandler(
    filename = function() paste0(input$microreact_project_name, "_metadata.csv"),
    content = function(file) {
      tryCatch({
        data <- create_microreact_data()
        writeLines(data$metadata_csv, file)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    }
  )

  output$download_microreact_tree <- downloadHandler(
    filename = function() paste0(input$microreact_project_name, "_tree.nwk"),
    content = function(file) {
      tryCatch({
        data <- create_microreact_data()
        writeLines(data$tree_newick, file)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    }
  )

  observeEvent(input$send_to_microreact, {
    # Local-only: create .microreact file and embed JSON in the local DB
    req(rv$tree_result)
    req(con())
    tryCatch({
      showNotification("Preparing Microreact data (local-only)...", type = "message", duration = 2)
      data <- create_microreact_data()
      timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
      safe_name <- gsub("[^A-Za-z0-9._-]", "_", input$microreact_project_name)
      filename <- paste0(safe_name, "_", timestamp, ".microreact")
      
      json_body <- build_microreact_project(data, input$microreact_project_name, input$microreact_description)

      # Save to local vendor viewer data dir (creates file for embedded viewer)
      viewer_data_dir <- file.path(app_dir, "viewer", "public", "data")
      dir.create(viewer_data_dir, recursive = TRUE, showWarnings = FALSE)
      out_path <- file.path(viewer_data_dir, filename)
      writeLines(json_body, out_path)

      # Embed in DB as an analysis result of type 'microreact'
      tryCatch({
        save_analysis_result(con(), analysis_name = input$microreact_project_name, analysis_type = "microreact", sample_ids = rv$selected_for_analysis, result_data = json_body, parameters = toJSON(list(saved_by = Sys.getenv("USERNAME"), created = Sys.time())))
      }, error = function(e) {
        warning("Failed to save microreact project to DB: ", e$message)
      })

      rv$current_microreact_url <- paste0("http://localhost:3000/?file=", URLencode(filename))
      rv$auto_open_local <- TRUE
      showNotification("Saved .microreact locally and embedded in DB (no server upload).", type = "message")
    }, error = function(e) {
      showNotification(paste("Error preparing Microreact data:", e$message), type = "error")
    })
  })

  # Render link and iframe
  output$microreact_link <- renderUI({
    if (!is.null(rv$current_microreact_url)) {
      tagList(p(strong("Project URL:")), tags$a(href = rv$current_microreact_url, target = "_blank", rv$current_microreact_url), br(), actionButton("open_microreact_iframe", "Load in Viewer Below", class = "btn-info btn-sm"))
    } else {
      p("No project created yet. Run analysis and click 'Send to Microreact'.")
    }
  })

  output$microreact_iframe <- renderUI({
    if (!is.null(rv$current_microreact_url) && ((!is.null(input$open_microreact_iframe) && input$open_microreact_iframe > 0) || isTRUE(rv$auto_open_local))) {
      tags$iframe(src = rv$current_microreact_url, width = "100%", height = "600px", frameborder = "0", style = "border: 1px solid #ddd; border-radius: 4px;")
    } else if (!is.null(rv$current_microreact_url)) {
      p("Click 'Load in Viewer Below' to display the visualization, or open the link in a new tab.")
    } else {
      p(em("Visualization will appear here after creating a project."))
    }
  })

  # Save local .microreact and set local viewer URL
  observeEvent(input$save_open_local, {
    req(rv$tree_result)
    tryCatch({
      data <- create_microreact_data()
      timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
      safe_name <- gsub("[^A-Za-z0-9._-]", "_", input$microreact_project_name)
      filename <- paste0(safe_name, "_", timestamp, ".microreact")
      
      json_body <- build_microreact_project(data, input$microreact_project_name, input$microreact_description)
      
      viewer_data_dir <- file.path(app_dir, "viewer", "public", "data")
      dir.create(viewer_data_dir, recursive = TRUE, showWarnings = FALSE)
      out_path <- file.path(viewer_data_dir, filename)
      writeLines(json_body, out_path)
      rv$current_microreact_url <- paste0("http://localhost:3000/?file=", URLencode(filename))
      rv$auto_open_local <- TRUE
      showNotification(paste("Saved .microreact to", out_path, "and set viewer URL."), type = "message")
    }, error = function(e) {
      showNotification(paste("Error saving local .microreact:", e$message), type = "error")
    })
  })

  output$microreact_history <- renderDT({
    if (nrow(rv$microreact_projects) > 0) {
      df <- rv$microreact_projects
      df$url <- paste0('<a href="', df$url, '" target="_blank">Open</a>')
      datatable(df, escape = FALSE, options = list(pageLength = 5), rownames = FALSE)
    } else {
      datatable(data.frame(Message = "No projects created in this session"), options = list(dom = 't'), rownames = FALSE)
    }
  })

  # ============================================================================
  # MICROREACT SERVER UPLOAD
  # ============================================================================

  # Check if team is configured (for conditional panel)
  output$has_team_configured <- reactive({
    !is.null(rv$microreact_team_id) && nchar(rv$microreact_team_id) > 0
  })
  outputOptions(output, "has_team_configured", suspendWhenHidden = FALSE)

  # Upload current analysis to Microreact server
  observeEvent(input$upload_to_microreact_server, {
    req(rv$tree_result)
    req(con())
    
    # Check for API token
    if (is.null(rv$microreact_token) || nchar(rv$microreact_token) == 0) {
      showNotification("Please configure your Microreact API token in the Settings tab first.", type = "error")
      return()
    }
    
    tryCatch({
      showNotification("Uploading to Microreact server...", type = "message", duration = 3)
      
      # Build the project JSON
      data <- create_microreact_data()
      json_body <- build_microreact_project(data, input$microreact_project_name, input$microreact_description)
      
      # DEBUG: Save the JSON being sent to a file for inspection
      debug_file <- file.path(app_dir, "debug_microreact_upload.json")
      writeLines(json_body, debug_file)
      message("DEBUG: Saved upload JSON to ", debug_file)
      
      # Upload to Microreact API
      result <- microreact_create_project(rv$microreact_token, json_body)
      
      if (result$success) {
        project_url <- result$data$url
        project_id <- result$data$id
        
        rv$current_microreact_url <- project_url
        rv$auto_open_local <- FALSE
        
        # Add to history
        new_row <- data.frame(
          name = input$microreact_project_name,
          url = project_url,
          id = project_id,
          created = as.character(Sys.time()),
          stringsAsFactors = FALSE
        )
        rv$microreact_projects <- rbind(rv$microreact_projects, new_row)
        
        # Share with team if configured and checkbox is checked
        if (isTRUE(input$share_with_team) && !is.null(rv$microreact_team_id) && nchar(rv$microreact_team_id) > 0) {
          share_result <- microreact_share_with_team(
            rv$microreact_token, 
            rv$microreact_team_id, 
            project_id, 
            input$team_share_role
          )
          if (share_result$success) {
            showNotification(paste("Project uploaded and shared with team as", input$team_share_role), type = "message")
          } else {
            showNotification(paste("Project uploaded but team sharing failed:", share_result$error), type = "warning")
          }
        } else {
          showNotification(paste("Project uploaded successfully! URL:", project_url), type = "message", duration = 10)
        }
        
        output$microreact_status <- renderPrint({
          cat("✓ Project uploaded to Microreact\n")
          cat("ID:", project_id, "\n")
          cat("URL:", project_url, "\n")
        })
        
      } else {
        showNotification(paste("Upload failed:", result$error), type = "error")
        output$microreact_status <- renderPrint({
          cat("✗ Upload failed\n")
          cat("Error:", result$error, "\n")
          if (!is.na(result$status)) {
            cat("HTTP Status:", result$status, "\n")
          }
        })
      }
      
    }, error = function(e) {
      showNotification(paste("Error uploading:", e$message), type = "error")
    })
  })

  # Upload existing .microreact file
  observeEvent(input$upload_existing_file, {
    req(input$microreact_file_upload)
    
    if (is.null(rv$microreact_token) || nchar(rv$microreact_token) == 0) {
      showNotification("Please configure your Microreact API token in the Settings tab first.", type = "error")
      return()
    }
    
    tryCatch({
      file_path <- input$microreact_file_upload$datapath
      showNotification("Uploading file to Microreact server...", type = "message", duration = 3)
      
      result <- microreact_create_project_from_file(rv$microreact_token, file_path)
      
      if (result$success) {
        project_url <- result$data$url
        project_id <- result$data$id
        
        rv$current_microreact_url <- project_url
        
        # Add to history
        new_row <- data.frame(
          name = input$microreact_file_upload$name,
          url = project_url,
          id = project_id,
          created = as.character(Sys.time()),
          stringsAsFactors = FALSE
        )
        rv$microreact_projects <- rbind(rv$microreact_projects, new_row)
        
        showNotification(paste("File uploaded successfully! URL:", project_url), type = "message", duration = 10)
        
        output$file_upload_status <- renderPrint({
          cat("✓ File uploaded to Microreact\n")
          cat("ID:", project_id, "\n")
          cat("URL:", project_url, "\n")
        })
        
      } else {
        showNotification(paste("Upload failed:", result$error), type = "error")
        output$file_upload_status <- renderPrint({
          cat("✗ Upload failed\n")
          cat("Error:", result$error, "\n")
        })
      }
      
    }, error = function(e) {
      showNotification(paste("Error uploading file:", e$message), type = "error")
    })
  })

  # ============================================================================
  # TEAM MANAGEMENT
  # ============================================================================

  # Save team ID
  observeEvent(input$save_team_id, {
    if (nchar(input$microreact_team_id) > 0) {
      rv$microreact_team_id <- input$microreact_team_id
      showNotification("Team ID saved for this session", type = "message")
      output$team_status <- renderPrint({
        cat("✓ Team ID saved:", input$microreact_team_id, "\n")
      })
    } else {
      showNotification("Please enter a team ID", type = "warning")
    }
  })

  # Create new team
  observeEvent(input$create_team, {
    req(input$new_team_name)
    
    if (is.null(rv$microreact_token) || nchar(rv$microreact_token) == 0) {
      showNotification("Please configure your Microreact API token first.", type = "error")
      return()
    }
    
    tryCatch({
      result <- microreact_create_team(rv$microreact_token, input$new_team_name)
      
      if (result$success) {
        team_id <- result$data$id
        rv$microreact_team_id <- team_id
        updateTextInput(session, "microreact_team_id", value = team_id)
        
        showNotification(paste("Team created! ID:", team_id), type = "message")
        output$create_team_status <- renderPrint({
          cat("✓ Team created successfully\n")
          cat("Team ID:", team_id, "\n")
          cat("Team Name:", input$new_team_name, "\n")
          cat("\nTeam ID has been saved automatically.\n")
        })
      } else {
        showNotification(paste("Failed to create team:", result$error), type = "error")
        output$create_team_status <- renderPrint({
          cat("✗ Failed to create team\n")
          cat("Error:", result$error, "\n")
        })
      }
    }, error = function(e) {
      showNotification(paste("Error creating team:", e$message), type = "error")
    })
  })

  # List team members
  observeEvent(input$list_team_members, {
    if (is.null(rv$microreact_token) || nchar(rv$microreact_token) == 0) {
      showNotification("Please configure your Microreact API token first.", type = "error")
      return()
    }
    
    if (is.null(rv$microreact_team_id) || nchar(rv$microreact_team_id) == 0) {
      # Try to use input directly
      if (nchar(input$microreact_team_id) > 0) {
        rv$microreact_team_id <- input$microreact_team_id
      } else {
        showNotification("Please save a team ID first.", type = "warning")
        return()
      }
    }
    
    tryCatch({
      result <- microreact_list_team_members(rv$microreact_token, rv$microreact_team_id)
      
      if (result$success) {
        output$team_members_list <- renderPrint({
          cat("Team Members:\n")
          cat("-------------\n")
          if (length(result$data) > 0) {
            members <- result$data
            if (is.data.frame(members)) {
              print(members)
            } else if (is.list(members)) {
              for (m in members) {
                if (is.list(m)) {
                  cat("- ", m$email %||% m$name %||% "Unknown", "\n")
                } else {
                  cat("- ", m, "\n")
                }
              }
            }
          } else {
            cat("No members found or empty response.\n")
          }
        })
        showNotification("Team members retrieved", type = "message")
      } else {
        output$team_members_list <- renderPrint({
          cat("✗ Failed to list members\n")
          cat("Error:", result$error, "\n")
        })
        showNotification(paste("Failed to list members:", result$error), type = "error")
      }
    }, error = function(e) {
      showNotification(paste("Error listing members:", e$message), type = "error")
    })
  })

  # Add team members
  observeEvent(input$add_team_members, {
    req(input$team_member_emails)
    
    if (is.null(rv$microreact_token) || nchar(rv$microreact_token) == 0) {
      showNotification("Please configure your Microreact API token first.", type = "error")
      return()
    }
    
    if (is.null(rv$microreact_team_id) || nchar(rv$microreact_team_id) == 0) {
      showNotification("Please save a team ID first.", type = "warning")
      return()
    }
    
    tryCatch({
      # Parse emails from textarea
      emails <- trimws(unlist(strsplit(input$team_member_emails, "\n")))
      emails <- emails[nchar(emails) > 0]
      
      if (length(emails) == 0) {
        showNotification("Please enter at least one email address", type = "warning")
        return()
      }
      
      result <- microreact_add_team_members(rv$microreact_token, rv$microreact_team_id, emails)
      
      if (result$success) {
        showNotification(paste("Added", length(emails), "member(s) to team"), type = "message")
        output$team_members_list <- renderPrint({
          cat("✓ Members added successfully:\n")
          for (e in emails) {
            cat("  -", e, "\n")
          }
        })
      } else {
        showNotification(paste("Failed to add members:", result$error), type = "error")
      }
    }, error = function(e) {
      showNotification(paste("Error adding members:", e$message), type = "error")
    })
  })

  # Remove team members
  observeEvent(input$remove_team_members, {
    req(input$team_remove_emails)
    
    if (is.null(rv$microreact_token) || nchar(rv$microreact_token) == 0) {
      showNotification("Please configure your Microreact API token first.", type = "error")
      return()
    }
    
    if (is.null(rv$microreact_team_id) || nchar(rv$microreact_team_id) == 0) {
      showNotification("Please save a team ID first.", type = "warning")
      return()
    }
    
    tryCatch({
      # Parse emails from textarea
      emails <- trimws(unlist(strsplit(input$team_remove_emails, "\n")))
      emails <- emails[nchar(emails) > 0]
      
      if (length(emails) == 0) {
        showNotification("Please enter at least one email address", type = "warning")
        return()
      }
      
      result <- microreact_remove_team_members(rv$microreact_token, rv$microreact_team_id, emails)
      
      if (result$success) {
        showNotification(paste("Removed", length(emails), "member(s) from team"), type = "message")
        output$team_members_list <- renderPrint({
          cat("✓ Members removed successfully:\n")
          for (e in emails) {
            cat("  -", e, "\n")
          }
        })
      } else {
        showNotification(paste("Failed to remove members:", result$error), type = "error")
      }
    }, error = function(e) {
      showNotification(paste("Error removing members:", e$message), type = "error")
    })
  })
}
