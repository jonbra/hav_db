ui_tab_settings <- function(){
  tabItem(tabName = "settings",
    fluidRow(
      box(title = "Microreact API Settings", status = "warning", solidHeader = TRUE, width = 6,
          p("Enter your Microreact API token to enable direct upload."),
          p("Get your token from: ", tags$a(href = "https://microreact.org/my-account/settings", target = "_blank", "https://microreact.org/my-account/settings")),
          hr(),
          passwordInput("microreact_api_token", "API Token:", value = ""),
          actionButton("save_api_token", "Save Token", class = "btn-primary"),
          actionButton("test_api_token", "Test Connection", class = "btn-info"),
          hr(),
          verbatimTextOutput("api_token_status")
      ),
      box(title = "About Microreact", status = "info", solidHeader = TRUE, width = 6,
          p("Microreact is a web application for visualization of phylogenetic trees with metadata."),
          p("Features:"),
          tags$ul(
            tags$li("Interactive phylogenetic tree visualization"),
            tags$li("Geographic mapping of samples"),
            tags$li("Timeline visualization"),
            tags$li("Metadata table with filtering")
          ),
          p("Learn more: ", tags$a(href = "https://microreact.org", target = "_blank", "https://microreact.org"))
      )
    )
  )
}
