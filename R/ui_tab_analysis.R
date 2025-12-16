ui_tab_analysis <- function(){
  tabItem(tabName = "analysis",
    fluidRow(
      box(title = "Selected Sequences for Analysis", status = "primary", solidHeader = TRUE, width = 12,
          DT::dataTableOutput("analysis_sequences"),
          actionButton("clear_analysis_selection", "Clear Selection", class = "btn-warning")
      )
    ),
    fluidRow(
      box(title = "Multiple Sequence Alignment", status = "info", solidHeader = TRUE, width = 6,
          selectInput("msa_method", "Alignment Method:", choices = c("Muscle", "ClustalW", "ClustalOmega")),
          actionButton("run_msa_btn", "Run Alignment", class = "btn-primary"),
          hr(),
          verbatimTextOutput("msa_stats"),
          downloadButton("download_alignment", "Download Alignment (FASTA)")
      ),
      box(title = "Phylogenetic Tree", status = "success", solidHeader = TRUE, width = 6,
          selectInput("tree_method", "Tree Method:", choices = c("Neighbor-Joining" = "nj", "UPGMA" = "upgma")),
          selectInput("dist_model", "Distance Model:", choices = c("K80", "K81", "F81", "F84", "T92", "TN93", "JC69", "raw")),
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
          numericInput("cluster_threshold", "Distance Threshold:", value = 0.03, min = 0, max = 1, step = 0.01),
          actionButton("run_cluster_btn", "Cluster Sequences", class = "btn-danger"),
          hr(),
          DT::dataTableOutput("cluster_results")
      )
    )
  )
}
