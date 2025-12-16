library(shiny)
library(shinydashboard)
library(DT)
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)

build_ui <- function(app_dir){
  # Source per-tab UI builders (requires `app_dir` to be defined by caller)
  source(file.path(app_dir, "R", "ui_tab_dashboard.R"))
  source(file.path(app_dir, "R", "ui_tab_browse.R"))
  source(file.path(app_dir, "R", "ui_tab_search.R"))
  source(file.path(app_dir, "R", "ui_tab_add.R"))
  source(file.path(app_dir, "R", "ui_tab_analysis.R"))
  source(file.path(app_dir, "R", "ui_tab_export.R"))
  source(file.path(app_dir, "R", "ui_tab_microreact.R"))
  source(file.path(app_dir, "R", "ui_tab_settings.R"))

  dashboardPage(
    dashboardHeader(title = "HAV Database"),

    dashboardSidebar(
      sidebarMenu(
        menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
        menuItem("Browse Data", tabName = "browse", icon = icon("database")),
        menuItem("Search", tabName = "search", icon = icon("search")),
        menuItem("Add Data", tabName = "add", icon = icon("plus")),
        menuItem("Analysis", tabName = "analysis", icon = icon("dna")),
        menuItem("Microreact", tabName = "microreact", icon = icon("project-diagram")),
        menuItem("Export", tabName = "export", icon = icon("download")),
        menuItem("Settings", tabName = "settings", icon = icon("cog"))
      )
    ),

    dashboardBody(
      tags$head(
        tags$style(HTML("\
          .sequence-display {\
            font-family: 'Courier New', monospace;\
            font-size: 12px;\
            word-wrap: break-word;\
            background-color: #f5f5f5;\
            padding: 10px;\
            border-radius: 4px;\
            max-height: 200px;\
            overflow-y: auto;\
          }\
          .info-box-icon { background-color: rgba(0,0,0,0.1) !important; }\
          .main-sidebar {\
            position: fixed !important;\
            top: 50px; /* adjust if header height differs */\
            left: 0;\
            height: calc(100vh - 50px);\
            overflow-y: auto;\
            z-index: 1000;\
          }\
          .content-wrapper, .main-footer {\
            margin-left: 230px !important; /* match sidebar width */\
          }\
          .main-header {\
            position: fixed !important;\
            top: 0;\
            left: 0;\
            right: 0;\
            z-index: 1100;\
            width: 100%;\
          }\
          .content-wrapper, .main-footer {\
            margin-top: 50px !important; /* same as header height */\
          }\
          .main-sidebar {\
            top: 50px; /* adjust if header height differs */\
          }\
        "))
      ),

      tabItems(
        ui_tab_dashboard(), ui_tab_browse(), ui_tab_search(), ui_tab_add(),
        ui_tab_analysis(), ui_tab_export(), ui_tab_microreact(), ui_tab_settings()
      )
    )
  )
}
