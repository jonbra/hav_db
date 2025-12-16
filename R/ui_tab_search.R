ui_tab_search <- function(){
  tabItem(tabName = "search",
    fluidRow(
      box(title = "Search Filters", status = "warning", solidHeader = TRUE, width = 4,
          selectInput("search_genotype", "Genotype:", choices = c("All" = ""), selected = ""),
          selectInput("search_variant", "Variant:", choices = c("All" = ""), selected = ""),
          textInput("search_patient", "Patient ID:", placeholder = "Search patient ID..."),
          selectInput("search_country", "Country:", choices = c("All" = ""), selected = ""),
          selectInput("search_location", "Location:", choices = c("All" = ""), selected = ""),
          numericInput("search_year_from", "Year From:", value = NA, min = 1990, max = 2100),
          numericInput("search_year_to", "Year To:", value = NA, min = 1990, max = 2100),
          numericInput("search_min_len", "Min Length:", value = NA, min = 0),
          numericInput("search_max_len", "Max Length:", value = NA, min = 0),
          actionButton("search_btn", "Search", class = "btn-primary", width = "100%"),
          hr(),
          actionButton("clear_search", "Clear Filters", width = "100%")
      ),
      checkboxGroupInput("search_filter_missing", "Hide samples missing:",
        choices = list(
          "Sampling date" = "sampling_date",
          "Country" = "geo_country",
          "Genotype" = "genotype"
        ),
        inline = FALSE
      ),
      box(title = "Search Results", status = "primary", solidHeader = TRUE, width = 8,
          DT::dataTableOutput("search_results"),
          hr(),
          fluidRow(
            column(4, actionButton("select_all_results", "Select All", class = "btn-info", width = "100%")),
            column(4, actionButton("deselect_all_results", "Deselect All", class = "btn-default", width = "100%")),
            column(4, actionButton("select_for_analysis", "Send to Analysis", class = "btn-success", width = "100%"))
          ),
          p(style = "margin-top: 10px;", textOutput("selection_count"))
      )
    )
  )
}
