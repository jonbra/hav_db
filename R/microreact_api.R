#' Microreact API Functions
#' 
#' Helper functions for interacting with the Microreact API.
#' Supports project creation, team management, and project sharing.
#' 
#' API Documentation: https://docs.microreact.org/api/creating-projects

library(httr)
library(jsonlite)

# Base URL for Microreact API
MICROREACT_API_BASE <- "https://microreact.org/api"

#' Validate API Token Format
#' 
#' @param token The API token to validate
#' @return TRUE if format looks valid, FALSE otherwise
microreact_validate_token_format <- function(token) {

if (is.null(token) || !is.character(token)) return(FALSE)
  # JWT tokens start with "eyJ"
  grepl("^eyJ", token) && nchar(token) > 50
}

#' Make API Request to Microreact
#' 
#' @param endpoint API endpoint (without base URL)
#' @param token API access token
#' @param body Request body (will be converted to JSON)
#' @param method HTTP method (POST, GET, etc.)
#' @return List with success status, response data, and any error message
microreact_api_request <- function(endpoint, token, body = NULL, method = "POST") {
  url <- paste0(MICROREACT_API_BASE, endpoint)
  
  headers <- add_headers(
    "Content-Type" = "application/json; charset=utf-8",
    "Access-Token" = token
  )
  
  tryCatch({
    if (method == "POST" && !is.null(body)) {
      # If body is already JSON string, use it directly
      if (is.character(body) && jsonlite::validate(body)) {
        json_body <- body
      } else {
        json_body <- toJSON(body, auto_unbox = TRUE)
      }
      response <- POST(url, headers, body = json_body, encode = "raw")
    } else if (method == "POST") {
      response <- POST(url, headers)
    } else {
      response <- GET(url, headers)
    }
    
    status <- status_code(response)
    content_text <- content(response, as = "text", encoding = "UTF-8")
    
    if (status >= 200 && status < 300) {
      result <- tryCatch(
        fromJSON(content_text),
        error = function(e) list(raw = content_text)
      )
      list(success = TRUE, data = result, status = status, error = NULL)
    } else {
      error_msg <- tryCatch(
        fromJSON(content_text)$message,
        error = function(e) content_text
      )
      list(success = FALSE, data = NULL, status = status, error = error_msg)
    }
  }, error = function(e) {
    list(success = FALSE, data = NULL, status = NA, error = e$message)
  })
}

# ============================================================================
# PROJECT MANAGEMENT
# ============================================================================

#' Create a Microreact Project from JSON
#' 
#' Upload a .microreact file (JSON) to create a new project.
#' 
#' @param token API access token
#' @param project_json JSON string or list containing the project definition
#' @return List with success status, project id, url, and any error
#' 
#' @examples
#' \dontrun{
#' result <- microreact_create_project(token, readLines("project.microreact"))
#' if (result$success) {
#'   cat("Project URL:", result$data$url, "\n")
#' }
#' }
microreact_create_project <- function(token, project_json) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  
  microreact_api_request("/projects/create", token, body = project_json)
}

#' Create Project from File
#' 
#' Read a .microreact file and upload it to create a project.
#' 
#' @param token API access token
#' @param file_path Path to the .microreact file
#' @return List with success status, project id, url, and any error
microreact_create_project_from_file <- function(token, file_path) {
  if (!file.exists(file_path)) {
    return(list(success = FALSE, data = NULL, error = "File not found"))
  }
  
  project_json <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  microreact_create_project(token, project_json)
}

#' Update an Existing Project
#' 
#' @param token API access token
#' @param project_id The project ID to update
#' @param project_json Updated project JSON
#' @return List with success status and any error
microreact_update_project <- function(token, project_id, project_json) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  
  endpoint <- paste0("/projects/", project_id, "/update")
  microreact_api_request(endpoint, token, body = project_json)
}

# ============================================================================
# TEAM MANAGEMENT (Experimental API)
# ============================================================================

#' Create a New Team
#' 
#' @param token API access token
#' @param team_name Name for the new team
#' @return List with success status, team id, and any error
microreact_create_team <- function(token, team_name) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  
  body <- list(name = team_name)
  microreact_api_request("/teams/create", token, body = body)
}

#' Add Members to a Team
#' 
#' @param token API access token
#' @param team_id The team ID (e.g., "53E42osSgysbaQGva5NciD")
#' @param emails Character vector of email addresses to add
#' @return List with success status and any error
microreact_add_team_members <- function(token, team_id, emails) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  
  if (!is.character(emails) || length(emails) == 0) {
    return(list(success = FALSE, data = NULL, error = "emails must be a non-empty character vector"))
  }
  
  body <- list(team = team_id, emails = as.list(emails))
  microreact_api_request("/teams/add-member", token, body = body)
}

#' List Team Members
#' 
#' @param token API access token
#' @param team_id The team ID
#' @return List with success status, member list, and any error
microreact_list_team_members <- function(token, team_id) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  
  body <- list(team = team_id)
  microreact_api_request("/teams/list-members", token, body = body)
}

#' Remove Members from a Team
#' 
#' @param token API access token
#' @param team_id The team ID
#' @param emails Character vector of email addresses to remove
#' @return List with success status and any error
microreact_remove_team_members <- function(token, team_id, emails) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  
  if (!is.character(emails) || length(emails) == 0) {
    return(list(success = FALSE, data = NULL, error = "emails must be a non-empty character vector"))
  }
  
  body <- list(team = team_id, emails = as.list(emails))
  microreact_api_request("/teams/remove-member", token, body = body)
}

# ============================================================================
# PROJECT SHARING
# ============================================================================

#' Share a Project with a Team
#' 
#' @param token API access token
#' @param team_id The team ID to share with
#' @param project_id The project ID to share
#' @param role Access role: "viewer", "editor", or "manager"
#' @return List with success status and any error
microreact_share_with_team <- function(token, team_id, project_id, role = "viewer") {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  
  valid_roles <- c("viewer", "editor", "manager")
  if (!role %in% valid_roles) {
    return(list(success = FALSE, data = NULL, error = paste("role must be one of:", paste(valid_roles, collapse = ", "))))
  }
  
  body <- list(team = team_id, project = project_id, role = role)
  microreact_api_request("/shares/add-team", token, body = body)
}

#' Unshare a Project from a Team
#' 
#' @param token API access token
#' @param team_id The team ID
#' @param project_id The project ID to unshare
#' @return List with success status and any error
microreact_unshare_from_team <- function(token, team_id, project_id) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  
  body <- list(team = team_id, project = project_id)
  microreact_api_request("/shares/remove-team", token, body = body)
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
  
#' Build Microreact Project JSON
#' 
#' Create a .microreact JSON structure from metadata and tree data.
#' 
#' @param metadata_csv CSV string of metadata
#' @param tree_newick Newick string of phylogenetic tree
#' @param project_name Name for the project
#' @param description Project description
#' @return JSON string suitable for upload
build_microreact_json <- function(metadata_csv, tree_newick, project_name = "Project", description = "") {
  timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  data_file_id <- paste0("data-", timestamp)
  tree_file_id <- paste0("tree-", timestamp)
  
  files_obj <- list()
  files_obj[[data_file_id]] <- list(
    name = "metadata.csv",
    format = "text/csv",
    blob = metadata_csv
  )
  files_obj[[tree_file_id]] <- list(
    name = "tree.nwk",
    format = "text/x-nh",
    blob = tree_newick
  )
  
  project <- list(
    meta = list(
      name = project_name,
      description = description
    ),
    files = files_obj,
    datasets = list(
      list(
        id = "dataset-1",
        file = data_file_id,
        idFieldName = "id"
      )
    ),
    trees = list(
      list(
        id = "tree-1",
        file = tree_file_id,
        labelField = "id"
      )
    )
  )
  
  toJSON(project, auto_unbox = TRUE)
}

#' Test API Connection
#' 
#' Attempt to verify the token works by making a minimal API call.
#' Note: There's no dedicated "ping" endpoint, so this creates a minimal
#' validation based on token format only.
#' 
#' @param token API access token
#' @return List with success status and message
microreact_test_connection <- function(token) {
  if (!microreact_validate_token_format(token)) {
    return(list(
      success = FALSE,
      message = "Token format is invalid. Microreact tokens are JWT format starting with 'eyJ'."
    ))
  }
  
  list(
    success = TRUE,
    message = "Token format is valid (JWT). Full validation occurs when creating a project."
  )
}
