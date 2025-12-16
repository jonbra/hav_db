library(DBI)
library(RSQLite)

ensure_provenance_table <- function(conn){
  DBI::dbExecute(conn, "CREATE TABLE IF NOT EXISTS provenance (\n+    id INTEGER PRIMARY KEY AUTOINCREMENT,\n+    action TEXT,\n+    actor TEXT,\n+    script TEXT,\n+    params TEXT,\n+    recorded_at TEXT DEFAULT (datetime('now'))\n+  )")
}

record_provenance <- function(conn, action, actor=NULL, script=NULL, params=NULL){
  ensure_provenance_table(conn)
  DBI::dbExecute(conn, "INSERT INTO provenance(action, actor, script, params) VALUES (?, ?, ?, ?)",
                params=list(action, actor, script, jsonlite::toJSON(params, auto_unbox=TRUE)))
}
