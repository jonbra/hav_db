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
          
          # Upload to Microreact Server (requires API token)
          h5(icon("cloud-upload-alt"), " Upload to Microreact.org"),
          actionButton("upload_to_microreact_server", "Upload to Microreact Server", 
                       class = "btn-success", icon = icon("upload"), width = "100%"),
          p(em("Requires API token configured in Settings tab"), style = "font-size: 11px; margin-top: 5px; color: #666;"),
          
          # Team sharing options (shown only if team is configured)
          conditionalPanel(
            condition = "output.has_team_configured",
            hr(),
            checkboxInput("share_with_team", "Share with team after upload", value = TRUE),
            selectInput("team_share_role", "Team role:", 
                        choices = c("Viewer" = "viewer", "Editor" = "editor", "Manager" = "manager"),
                        selected = "viewer")
          ),
          
          hr(),
          h5(icon("download"), " Manual Upload Options:"),
          downloadButton("download_microreact_file", "Download Files (ZIP)", class = "btn-info", style = "width: 100%;"),
          p(em("Downloads metadata.csv + tree.nwk - upload both to microreact.org/upload"), style = "font-size: 11px; margin-top: 5px;"),
          fluidRow(
            column(6, downloadButton("download_microreact_csv", "Download CSV", class = "btn-default", style = "width: 100%;")),
            column(6, downloadButton("download_microreact_tree", "Download Tree", class = "btn-default", style = "width: 100%;"))
          ),
          hr(),
          h5(icon("desktop"), " Local Options:"),
          actionButton("send_to_microreact", "Save Locally", class = "btn-default", icon = icon("save"), width = "100%"),
          p(em("Saves .microreact file locally and to database"), style = "font-size:11px; margin-top:5px;"),
          actionButton("save_open_local", "Save & Open in Local Viewer", class = "btn-warning", width = "100%"),
          p(em("Writes .microreact into viewer/public/data and opens embedded viewer"), style = "font-size:11px; margin-top:5px;"),
          verbatimTextOutput("microreact_status")
      ),
      box(title = "Microreact Visualization", status = "success", solidHeader = TRUE, width = 6,
          p("After uploading to Microreact, the visualization will appear below."),
          uiOutput("microreact_link"),
          hr(),
          uiOutput("microreact_iframe")
      )
    ),
    fluidRow(
      box(title = "Upload .microreact File", status = "warning", solidHeader = TRUE, width = 6,
          p("Upload an existing .microreact file to microreact.org"),
          fileInput("microreact_file_upload", "Choose .microreact file:", 
                    accept = c(".microreact", ".json")),
          actionButton("upload_existing_file", "Upload File to Server", 
                       class = "btn-primary", icon = icon("cloud-upload-alt")),
          verbatimTextOutput("file_upload_status")
      ),
      box(title = "Previous Microreact Projects", status = "info", solidHeader = TRUE, width = 6,
          p("Your recent Microreact projects from this session:"),
          DT::dataTableOutput("microreact_history")
      )
    )
  )
}
