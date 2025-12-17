library(shiny)

ui <- fluidPage(
  titlePanel("Shiny deployment test"),
  mainPanel(
    h2("Shiny test"),
    p("Current time:"),
    textOutput("time")
  )
)

server <- function(input, output, session) {
  output$time <- renderText(Sys.time())
}

# Host/port from environment with sensible defaults
host <- Sys.getenv("SHINY_HOST", "0.0.0.0")
port <- as.integer(Sys.getenv("SHINY_PORT", "3838"))

shiny::runApp(list(ui = ui, server = server), host = host, port = port)
