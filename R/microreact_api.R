#' Microreact API Functions
#' 
#' Helper functions for interacting with the Microreact API.
#' Supports project creation, team management, and project sharing.
#' 
#' API Documentation: https://docs.microreact.org/api/creating-projects

library(httr)
library(jsonlite)
library(base64enc)

# Base URL for Microreact API
MICROREACT_API_BASE <- "https://microreact.org/api"

#' Validate API Token Format
#' 
#' @param token The API token to validate
#' @return TRUE if format looks valid, FALSE otherwise
microreact_validate_token_format <- function(token) {
  if (is.null(token) || !is.character(token)) return(FALSE)
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

microreact_create_project <- function(token, project_json) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  microreact_api_request("/projects/create", token, body = project_json)
}

microreact_create_project_from_file <- function(token, file_path) {
  if (!file.exists(file_path)) {
    return(list(success = FALSE, data = NULL, error = "File not found"))
  }
  project_json <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  microreact_create_project(token, project_json)
}

microreact_update_project <- function(token, project_id, project_json) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  endpoint <- paste0("/projects/", project_id, "/update")
  microreact_api_request(endpoint, token, body = project_json)
}

# ============================================================================
# TEAM MANAGEMENT
# ============================================================================

microreact_create_team <- function(token, team_name) {
  if (!microreact_validate_token_format(token)) {
    return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  }
  body <- list(name = team_name)
  microreact_api_request("/teams/create", token, body = body)
}

microreact_add_team_members <- function(token, team_id, emails) {
  if (!microreact_validate_token_format(token)) return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  if (!is.character(emails) || length(emails) == 0) return(list(success = FALSE, data = NULL, error = "emails must be a non-empty character vector"))
  body <- list(team = team_id, emails = as.list(emails))
  microreact_api_request("/teams/add-member", token, body = body)
}

microreact_list_team_members <- function(token, team_id) {
  if (!microreact_validate_token_format(token)) return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  body <- list(team = team_id)
  microreact_api_request("/teams/list-members", token, body = body)
}

microreact_remove_team_members <- function(token, team_id, emails) {
  if (!microreact_validate_token_format(token)) return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  if (!is.character(emails) || length(emails) == 0) return(list(success = FALSE, data = NULL, error = "emails must be a non-empty character vector"))
  body <- list(team = team_id, emails = as.list(emails))
  microreact_api_request("/teams/remove-member", token, body = body)
}

# ============================================================================
# PROJECT SHARING
# ============================================================================

microreact_share_with_team <- function(token, team_id, project_id, role = "viewer") {
  if (!microreact_validate_token_format(token)) return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  valid_roles <- c("viewer", "editor", "manager")
  if (!role %in% valid_roles) return(list(success = FALSE, data = NULL, error = paste("role must be one of:", paste(valid_roles, collapse = ", "))))
  body <- list(team = team_id, project = project_id, role = role)
  microreact_api_request("/shares/add-team", token, body = body)
}

microreact_unshare_from_team <- function(token, team_id, project_id) {
  if (!microreact_validate_token_format(token)) return(list(success = FALSE, data = NULL, error = "Invalid token format"))
  body <- list(team = team_id, project = project_id)
  microreact_api_request("/shares/remove-team", token, body = body)
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

build_microreact_payload_from_paths <- function(metadata_csv_path, tree_nwk_path,
                                                project_name = "HAV_Analysis",
                                                description = "HAV phylogenetic analysis",
                                                include_image_path = NULL) {
  if (!file.exists(metadata_csv_path)) stop("metadata CSV not found: ", metadata_csv_path)
  if (!file.exists(tree_nwk_path)) stop("tree file not found: ", tree_nwk_path)

  ts <- format(Sys.time(), "%Y%m%d%H%M%S")
  data_id <- paste0("data", ts)
  tree_id <- paste0("tree", ts)

  csv_b64 <- base64enc::base64encode(metadata_csv_path)
  tree_b64 <- base64enc::base64encode(tree_nwk_path)

  size_csv <- as.integer(file.info(metadata_csv_path)$size)
  size_tree <- as.integer(file.info(tree_nwk_path)$size)

  files <- list()
  files[[data_id]] <- list(
    id = data_id,
    name = basename(metadata_csv_path),
    format = "text/csv",
    type = "data",
    size = size_csv,
    blob = paste0("data:text/csv;base64,", csv_b64)
  )
  files[[tree_id]] <- list(
    id = tree_id,
    name = basename(tree_nwk_path),
    format = "text/x-nh",
    type = "tree",
    size = size_tree,
    blob = paste0("data:application/octet-stream;base64,", tree_b64)
  )

  meta_time_iso <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")

  payload <- list(
    schema = "https://microreact.org/schema/v1.json",
    meta = list(
      name = project_name,
      description = description,
      timestamp = meta_time_iso
    ),
    files = files,
    datasets = list(
      `dataset-1` = list(
        id = "dataset-1",
        file = data_id,
        idFieldName = "id"
      )
    ),
    trees = list(
      `tree-1` = list(
        id = "tree-1",
        title = "Tree",
        file = tree_id,
        labelsField = "id",
        type = "rc",
        alignLabels = TRUE,
        showLabels = TRUE,
        showLeafLabels = TRUE,
        showShapes = TRUE,
        showShapeBorders = TRUE,
        showInternalLabels = TRUE,
        showPiecharts = TRUE,
        showEdges = TRUE,
        nodeSize = 14,
        fontSize = 16,
        controls = TRUE,
        blocks = list(),
        showBlockHeaders = FALSE,
        blockSize = 14
      )
    ),
    styles = list(
      nodes = list(
        fill = "#1f77b4",
        stroke = "#000000",
        strokeWidth = 1
      ),
      edges = list(
        stroke = "#999999",
        strokeWidth = 1
      ),
      labels = list(
        font = "Arial",
        size = 14,
        color = "#333333"
      )
    )
  )

  if (!is.null(include_image_path) && file.exists(include_image_path)) {
    img_b64 <- base64enc::base64encode(include_image_path)
    payload$meta$image <- paste0("data:image/png;base64,", img_b64)
  }

  toJSON(payload, auto_unbox = TRUE, pretty = TRUE)
}

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
