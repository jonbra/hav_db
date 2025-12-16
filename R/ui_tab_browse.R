ui_tab_browse <- function(){
  tabItem(tabName = "browse",
    fluidRow(
      box(title = "All Sequences", status = "primary", solidHeader = TRUE, width = 12,
          DT::dataTableOutput("all_sequences_table"),
          hr(),
          h4("Selected Sequence Details"),
          verbatimTextOutput("selected_sequence_info"),
          div(class = "sequence-display", verbatimTextOutput("selected_sequence_seq"))
      )
    )
  )
}
