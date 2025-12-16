register_microreact_server <- function(input, output, session, rv, con, app_dir){
  # Microreact Settings handlers
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

  observeEvent(input$test_api_token, {
    req(input$microreact_api_token)
    token <- input$microreact_api_token
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

  # Helper to build microreact data
  create_microreact_data <- function(){
    req(rv$tree_result)
    req(con())
    sample_ids <- rv$selected_for_analysis
    metadata <- get_sequences_with_metadata(con(), sample_ids)
    tree_tip_labels <- rv$tree_result$tip.label
    id_to_tree_label <- function(id) { gsub(" ", "_", id) }
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
    tree_newick <- write.tree(rv$tree_result)
    temp_csv <- tempfile(fileext = ".csv")
    write.csv(microreact_meta, temp_csv, row.names = FALSE, na = "")
    meta_csv <- paste(readLines(temp_csv, warn = FALSE), collapse = "\n")
    unlink(temp_csv)
    list(metadata_csv = meta_csv, tree_newick = tree_newick, metadata_df = microreact_meta)
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
      data_file_id <- paste0("data-", timestamp)
      tree_file_id <- paste0("tree-", timestamp)
      files_obj <- list()
      files_obj[[data_file_id]] <- list(name = "metadata.csv", format = "text/csv", blob = data$metadata_csv)
      files_obj[[tree_file_id]] <- list(name = "tree.nwk", format = "text/x-nh", blob = data$tree_newick)
      microreact_project <- list(
        meta = list(name = input$microreact_project_name, description = input$microreact_description),
        files = files_obj,
        datasets = list(list(id = "dataset-1", file = data_file_id, idFieldName = "id")),
        trees = list(list(id = "tree-1", file = tree_file_id, labelField = "id")),
        initial_view = list(labelField = "country", colourField = "country", colourColumn = "country__colour", showMap = TRUE)
      )
      json_body <- toJSON(microreact_project, auto_unbox = TRUE)

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
      data_file_id <- paste0("data-", timestamp)
      tree_file_id <- paste0("tree-", timestamp)
      files_obj <- list()
      files_obj[[data_file_id]] <- list(name = "metadata.csv", format = "text/csv", blob = data$metadata_csv)
      files_obj[[tree_file_id]] <- list(name = "tree.nwk", format = "text/x-nh", blob = data$tree_newick)
      microreact_project <- list(meta = list(name = input$microreact_project_name, description = input$microreact_description), files = files_obj, datasets = list(list(id = "dataset-1", file = data_file_id, idFieldName = "id")), trees = list(list(id = "tree-1", file = tree_file_id, labelField = "id")), initial_view = list(labelField = "country", colourField = "country", colourColumn = "country__colour", showMap = TRUE))
      json_body <- toJSON(microreact_project, auto_unbox = TRUE)
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
}
