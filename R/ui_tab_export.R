ui_tab_export <- function(){
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
  )
}
