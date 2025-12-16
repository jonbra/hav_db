library(DBI)
library(RSQLite)

# QC helper functions for the database
connect_db <- function(db_path){
  DBI::dbConnect(RSQLite::SQLite(), db_path)
}

check_missing_sample_fields <- function(conn){
  df <- DBI::dbGetQuery(conn, "SELECT id, sample_id FROM samples WHERE sample_id IS NULL OR sample_id = ''")
  if(nrow(df) == 0) return(list(ok=TRUE, msg='no missing sample_id'))
  list(ok=FALSE, msg=paste0(nrow(df),' samples with missing sample_id'))
}

check_duplicate_sequence_ids <- function(conn){
  df <- DBI::dbGetQuery(conn, "SELECT sequence_id, COUNT(*) AS cnt FROM sequences GROUP BY sequence_id HAVING cnt>1")
  if(nrow(df) == 0) return(list(ok=TRUE, msg='no duplicate sequence_id'))
  list(ok=FALSE, msg=paste0(nrow(df),' duplicate sequence_id(s)'), details=df)
}

check_date_ranges <- function(conn){
  df <- tryCatch(DBI::dbGetQuery(conn, "SELECT id, collection_date FROM samples WHERE collection_date IS NOT NULL"), error=function(e) NULL)
  if(is.null(df) || nrow(df)==0) return(list(ok=TRUE, msg='no collection dates to check'))
  # basic parse attempt
  bad <- df[is.na(as.Date(df$collection_date, format='%Y-%m-%d')), ]
  if(nrow(bad)==0) return(list(ok=TRUE, msg='dates parse successfully'))
  list(ok=FALSE, msg=paste0(nrow(bad),' samples with unparseable dates'), details=bad)
}

run_all_checks <- function(db_path, verbose=TRUE){
  conn <- connect_db(db_path)
  on.exit(DBI::dbDisconnect(conn))

  res <- list()
  res$missing_sample_fields <- check_missing_sample_fields(conn)
  res$duplicate_seq_ids <- check_duplicate_sequence_ids(conn)
  res$date_checks <- check_date_ranges(conn)

  if(verbose){
    for(n in names(res)){
      cat(sprintf('%s: %s\n', n, res[[n]]$msg))
    }
  }
  res
}
