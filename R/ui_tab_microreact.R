ui_tab_microreact <- function(){
  tabItem(tabName = "microreact",
    fluidRow(
      box(title = "Send to Microreact", status = "primary", solidHeader = TRUE, width = 6,
          p("Create a Microreact visualization from your phylogenetic tree and metadata."),
          p(strong("Requirements:"), "Run alignment and build a tree first in the Analysis tab."),
          hr(),
          textInput("microreact_project_name", "Project Name:", value = paste0("HAV_Analysis_", format(Sys.Date(), "%Y%m%d"))),
          textAreaInput("microreact_description", "Description:", rows = 2, value = "HAV phylogenetic analysis"),
          hr(),
          actionButton("send_to_microreact", "Send to Microreact", class = "btn-success", icon = icon("upload"), width = "100%"),
          br(), br(),
          h5("Manual Upload Options:"),
          downloadButton("download_microreact_file", "Download Files (ZIP)", class = "btn-info", style = "width: 100%;"),
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
  )
}
