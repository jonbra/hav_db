library(DBI)
library(RSQLite)

con <- dbConnect(SQLite(), "data/hav.db")

# Check for QCMD and test samples
res <- dbGetQuery(con, "SELECT sample_id, geo_country, geo_location FROM metadata 
                         WHERE sample_id LIKE '%QCMD%' 
                         OR sample_id LIKE '%596%' 
                         OR sample_id LIKE '%953%'
                         OR sample_id LIKE '%Gt %'")
print("Samples matching QCMD/596/953/Gt:")
print(res)

dbDisconnect(con)
