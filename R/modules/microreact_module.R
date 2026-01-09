register_microreact_server <- function(input, output, session, rv, con, app_dir){
  # Source the API functions
  source(file.path(app_dir, "R", "microreact_api.R"), local = TRUE)
  
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
      timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
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
