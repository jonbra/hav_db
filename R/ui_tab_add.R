ui_tab_add <- function(){
  tabItem(tabName = "add",
    fluidRow(
      box(title = "Add Single Sequence", status = "success", solidHeader = TRUE, width = 6,
          textInput("new_sample_id", "Sample ID:", placeholder = "e.g., SAMPLE_001"),
          textAreaInput("new_sequence", "Sequence:", rows = 5, placeholder = "Paste DNA sequence here..."),
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
          p("CSV must contain 'sample_id' column. Other columns: sampling_date, sample_year, genotype, \n                 variant, patient_id, geo_location, geo_country, source, transmission_route, comment"),
          hr(),
          actionButton("import_btn", "Import", class = "btn-primary", width = "100%"),
          verbatimTextOutput("import_result")
      )
    )
  )
}
