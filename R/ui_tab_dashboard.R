ui_tab_dashboard <- function(){
  tabItem(tabName = "dashboard",
    fluidRow(
      infoBoxOutput("total_sequences", width = 4),
      infoBoxOutput("total_viruses", width = 4),
      infoBoxOutput("date_range", width = 4)
    ),
    fluidRow(
      box(title = "Sequences by Virus", status = "primary", solidHeader = TRUE,
          plotlyOutput("virus_plot", height = 300), width = 6),
      box(title = "Sequences by Location", status = "primary", solidHeader = TRUE,
          plotlyOutput("location_plot", height = 300), width = 6)
    ),
    fluidRow(
      box(title = "Recent Sequences", status = "info", solidHeader = TRUE,
          DT::dataTableOutput("recent_sequences"), width = 12)
    )
  )
}
