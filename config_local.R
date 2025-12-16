## Local application configuration (safe defaults)
# This file is sourced by app.R. It prefers environment overrides and
# falls back to sensible defaults when `app_dir` is not available.

`%||%` <- function(x, y) if (is.null(x)) y else x

DB_PATH <- Sys.getenv("HAV_DB_PATH", unset = "")
if (identical(DB_PATH, "") || is.null(DB_PATH)) {
  if (exists("app_dir") && !is.null(app_dir) && nzchar(app_dir)) {
    DB_PATH <- file.path(app_dir, "data", "hav.db")
  } else {
    DB_PATH <- file.path(dirname(sys.frame(1)$ofile %||% "."), "data", "hav.db")
  }
}

# Directory used to store local viewer data (.microreact files)
VIEWER_DIR <- Sys.getenv("HAV_VIEWER_DIR", unset = file.path(dirname(sys.frame(1)$ofile %||% "."), "viewer", "vendor", "public"))
