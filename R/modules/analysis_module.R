register_analysis_server <- function(input, output, session, rv, con){
  # Renders the selected sequences table
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
    cat("Conserved positions:", stats$conserved_positions, "(", stats$conservation_pct, "% )\n")
    cat("Gap positions:", stats$gap_positions, "\n")
  })

  output$download_alignment <- downloadHandler(
    filename = function() paste0("alignment_", Sys.Date(), ".fasta"),
    content = function(file) {
      req(rv$msa_result)
      export_alignment_fasta(rv$msa_result, file)
    }
  )

  # BLAST-like nearest-match search (R-based)
observeEvent(input$find_closest_btn, {
  req(con())
  # read query
  qset <- NULL
  if (!is.null(input$blast_query_file) && nzchar(input$blast_query_file$datapath)) {
    qfile <- input$blast_query_file$datapath
    qset <- tryCatch(readDNAStringSet(qfile, format = "fasta"), error = function(e) NULL)
  } else if (!is.null(input$blast_query_seq) && nzchar(input$blast_query_seq)) {
    txt <- input$blast_query_seq
    if (grepl("^>", trimws(txt))) {
      tf <- tempfile(fileext = ".fa")
      writeLines(txt, tf)
      qset <- tryCatch(readDNAStringSet(tf, format = "fasta"), error = function(e) NULL)
    } else {
      qset <- DNAStringSet(toupper(gsub("\\s+", "", txt)))
      names(qset) <- "query_1"
    }
  } else {
    showNotification("Provide a query sequence or upload a FASTA file", type = "warning")
    return()
  }

  # Run the R-based nearest-match function
  tryCatch({
    max_hits <- as.integer(input$blast_max_hits %||% 10)
    min_id <- as.numeric(input$blast_min_identity %||% 0)
    res <- run_blast(query = qset, con = con(), max_hits = max_hits, min_identity = min_id)
    if (is.null(res) || is.null(res$hits) || nrow(res$hits) == 0) {
      output$blast_hits <- renderDT({ datatable(data.frame(Message = "No hits found"), options = list(dom = 't')) })
      rv$blast_hits <- NULL
      showNotification("No hits found", type = "message")
      return()
    }
    rv$blast_hits <- res
    output$blast_hits <- renderDT({
      datatable(res$hits, options = list(pageLength = 25), rownames = FALSE)
    })
    showNotification("Nearest-match search complete", type = "message")
  }, error = function(e) {
    showNotification(paste("Search failed:", e$message), type = "error")
  })
})

output$download_blast_results <- downloadHandler(
  filename = function() paste0("blast_hits_", Sys.Date(), ".fasta"),
  content = function(file) {
    req(rv$blast_hits)
    if (!is.null(rv$blast_hits$sequences) && length(rv$blast_hits$sequences) > 0) {
      Biostrings::writeXStringSet(rv$blast_hits$sequences, filepath = file, format = "fasta")
    } else {
      writeLines(">no_sequences", file)
    }
  }
)

  # Tree
  observeEvent(input$run_tree_btn, {
    req(rv$msa_result)

    withProgress(message = "Building tree...", {
      tryCatch({
        if (input$tree_method == "nj") {
          rv$tree_result <- build_nj_tree(rv$msa_result, input$dist_model)
        } else if (input$tree_method == "upgma") {
          rv$tree_result <- build_upgma_tree(rv$msa_result, input$dist_model)
        } else if (input$tree_method == "iqtree") {
          # Export alignment to fasta and run IQ-TREE wrapper
          fa <- tempfile(fileext = ".fa")
          export_alignment_fasta(rv$msa_result, fa)
          pref <- tempfile("iqtree")
          rv$tree_result <- run_iqtree(fa, prefix = pref, threads = as.integer(input$tree_threads %||% 1))
        } else {
          stop("Unsupported tree method")
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

  # Analysis persistence (save / load / delete)
  analysis_dir <- file.path(dirname(sys.frame(1)$ofile %||% "."), "analysis")
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
    infos2 <- Filter(Negate(is.null), infos)
    if (length(infos2) == 0) return(data.frame())
    do.call(rbind, infos2)
  }

  output$stored_analyses <- renderDT({
    df <- list_saved_analyses()
    if (is.null(df) || nrow(df) == 0) return(datatable(data.frame(Message = "No saved analyses"), options = list(dom = 't')))
    disp <- df[, c("file", "name", "type", "created", "description")]
    datatable(disp, selection = list(mode = 'single', target = 'row'), rownames = FALSE, escape = TRUE, options = list(pageLength = 10))
  })

  observeEvent(input$save_analysis, {
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
    output$stored_analyses <- renderDT({
      df <- list_saved_analyses()
      if (is.null(df) || nrow(df) == 0) return(datatable(data.frame(Message = "No saved analyses"), options = list(dom = 't')))
      disp <- df[, c("file", "name", "type", "created", "description")]
      datatable(disp, selection = list(mode = 'single', target = 'row'), rownames = FALSE, escape = TRUE, options = list(pageLength = 10))
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
    loaded_tree <- obj$tree
    if (!is.null(loaded_tree)) {
      if (is.character(loaded_tree)) {
        loaded_tree2 <- tryCatch(ape::read.tree(text = loaded_tree), error = function(e) NULL)
      } else if (inherits(loaded_tree, "phylo")) {
        loaded_tree2 <- loaded_tree
      } else if (is.list(loaded_tree) && !is.null(loaded_tree$tip.label)) {
        class(loaded_tree) <- c("phylo", class(loaded_tree))
        loaded_tree2 <- loaded_tree
      } else {
        loaded_tree2 <- NULL
      }
    } else {
      loaded_tree2 <- NULL
    }
    rv$msa_result <- obj$msa
    rv$tree_result <- loaded_tree2
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
    output$stored_analyses <- renderDT({
      df <- list_saved_analyses()
      if (is.null(df) || nrow(df) == 0) return(datatable(data.frame(Message = "No saved analyses"), options = list(dom = 't')))
      disp <- df[, c("file", "name", "type", "created", "description")]
      datatable(disp, selection = list(mode = 'single', target = 'row'), rownames = FALSE, escape = TRUE, options = list(pageLength = 10))
    })
  })
}
