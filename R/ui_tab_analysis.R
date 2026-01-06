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
          selectInput("msa_method", "Alignment Method:", choices = c("MAFFT" = "mafft"), selected = "mafft"),
          actionButton("run_msa_btn", "Run Alignment", class = "btn-primary"),
          verbatimTextOutput("msa_stats"),
          downloadButton("download_alignment", "Download Alignment (FASTA)"),
          hr(),
          textAreaInput("blast_query_seq", "Query sequence (FASTA or raw sequence):", rows = 4, placeholder = "Paste sequence or FASTA here"),
          fileInput("blast_query_file", "Or upload FASTA file:", accept = c(".fa", ".fasta")),
          numericInput("blast_max_hits", "Max hits per query:", value = 10, min = 1, step = 1),
          numericInput("blast_min_identity", "Min percent identity:", value = 70, min = 0, max = 100, step = 1),
          actionButton("find_closest_btn", "Find Closest Sequences", class = "btn-primary"),
          downloadButton("download_blast_results", "Download BLAST Hits (FASTA)"),
          DT::dataTableOutput("blast_hits")
      ),
      box(title = "Phylogenetic Tree", status = "success", solidHeader = TRUE, width = 6,
          selectInput("tree_method", "Tree Method:", choices = c("Neighbor-Joining" = "nj", "UPGMA" = "upgma", "IQ-TREE (external)" = "iqtree")),
          # Distance model only for NJ/UPGMA
          conditionalPanel(
            condition = "input.tree_method != 'iqtree'",
            selectInput("dist_model", "Distance Model:", choices = c("K80", "K81", "F81", "F84", "T92", "TN93", "JC69", "raw"))
          ),
          # IQ-TREE specific options
          conditionalPanel(
            condition = "input.tree_method == 'iqtree'",
            textInput("iqtree_model", "Substitution Model:", value = "GTR+G+I", 
                      placeholder = "e.g., GTR+G+I, HKY+G, TIM+I, MFP (auto)"),
            numericInput("iqtree_bootstrap", "Ultrafast Bootstrap Replicates (-B):", 
                         value = 1000, min = 0, max = 100000, step = 100),
            helpText("Set to 0 to disable bootstrap analysis."),
            numericInput("iqtree_threads", "CPU Threads (-T):", value = 2, min = 1, max = 64, step = 1)
          ),
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
