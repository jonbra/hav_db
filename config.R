# =============================================================================
# HAV Database Configuration
# =============================================================================
# This file contains configuration settings for the HAV database.
# When moving to a shared drive, update the DB_PATH variable below.

# Database path - UPDATE THIS when moving to shared drive
# Example for shared drive: "N:/shared/hav_db/data/hav.db"
DB_PATH <- file.path(dirname(sys.frame(1)$ofile %||% "."), "data", "hav.db")

# Alternative: Set an explicit path (uncomment and modify when needed)
# DB_PATH <- "C:/Users/jonr/OneDrive - Folkehelseinstituttet/Prosjekter/hav_db/data/hav.db"

# Null coalescing operator if not available
`%||%` <- function(x, y) if (is.null(x)) y else x
