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
    ),
    fluidRow(
      box(title = "Team Management (Experimental)", status = "primary", solidHeader = TRUE, width = 6,
          p("Configure team sharing for Microreact projects."),
          p(em("Note: Team features require an API token to be configured above.")),
          hr(),
          h5(icon("users"), " Current Team"),
          textInput("microreact_team_id", "Team ID:", value = "", placeholder = "e.g., 53E42osSgysbaQGva5NciD"),
          actionButton("save_team_id", "Save Team ID", class = "btn-primary"),
          actionButton("list_team_members", "List Members", class = "btn-info"),
          verbatimTextOutput("team_status"),
          hr(),
          h5(icon("user-plus"), " Create New Team"),
          textInput("new_team_name", "Team Name:", value = "", placeholder = "Enter team name"),
          actionButton("create_team", "Create Team", class = "btn-success"),
          verbatimTextOutput("create_team_status")
      ),
      box(title = "Team Members", status = "info", solidHeader = TRUE, width = 6,
          h5(icon("user-plus"), " Add Members"),
          textAreaInput("team_member_emails", "Email addresses (one per line):", rows = 3, 
                        placeholder = "user1@example.com\nuser2@example.com"),
          actionButton("add_team_members", "Add Members", class = "btn-success"),
          hr(),
          h5(icon("user-minus"), " Remove Members"),
          textAreaInput("team_remove_emails", "Email addresses to remove (one per line):", rows = 3),
          actionButton("remove_team_members", "Remove Members", class = "btn-danger"),
          hr(),
          verbatimTextOutput("team_members_list")
      )
    )
  )
}
