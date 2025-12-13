# HAV Database

A lightweight local database for storing and analyzing HAV sequences with metadata.

## Features

- **SQLite database** - Fast, portable, single-file storage
- **Shiny web interface** - Browse, search, and analyze sequences
- **R integration** - Full command-line access via R functions
- **Analysis tools** - Multiple sequence alignment, phylogenetic trees, clustering
- **Import/Export** - FASTA and CSV support

## Quick Start

### 1. Install Required Packages

```r
# CRAN packages
install.packages(c(
  "DBI", "RSQLite", "dbplyr", "dplyr",
  "seqinr", "ape", "jsonlite",
  "shiny", "shinydashboard", "DT", "plotly", "ggplot2"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "msa"))
```

### 2. Initialize the Database

```r
setwd("path/to/hav_db")
source("setup_database.R")
```

### 3. Launch the Shiny App

```r
shiny::runApp("path/to/hav_db")
```

Or from within RStudio, open `app.R` and click "Run App".

## Directory Structure

```
hav_db/
├── app.R                    # Shiny application
├── config.R                 # Configuration (database path)
├── setup_database.R         # Database initialization script
├── README.md
├── R/
│   ├── db_functions.R       # Database CRUD operations
│   ├── sequence_functions.R # FASTA import/export
│   └── analysis_functions.R # MSA, phylogeny, clustering
├── data/
│   └── hav.db               # SQLite database (created by setup)
└── www/                     # Static assets for Shiny
```

## Command-Line Usage

### Load Functions

```r
source("R/db_functions.R")
source("R/sequence_functions.R")
source("R/analysis_functions.R")

# Connect to database
con <- get_db_connection("data/hav.db")
```

### Add Sequences

```r
# Single sequence with metadata
add_sequence_with_metadata(con,
  sample_id = "SAMPLE_001",
  sequence = "ATGCGATCGATCGATCG...",
  virus_name = "Hepatitis A Virus",
  sampling_date = "2025-01-15",
  geo_location = "Oslo"
)

# Import from FASTA
import_fasta_to_db(con, "sequences.fasta")

# Import with metadata CSV
metadata <- read.csv("metadata.csv")
import_fasta_to_db(con, "sequences.fasta", metadata_df = metadata)
```

### Search and Query

```r
# Get all sequences
all_seqs <- get_all_sequences(con)

# Search with filters
results <- search_sequences(con,
  virus_name = "Hepatitis A",
  date_from = "2024-01-01",
  date_to = "2024-12-31"
)

# Get sequences with metadata
data <- get_sequences_with_metadata(con, c("SAMPLE_001", "SAMPLE_002"))
```

### Run Analysis

```r
# Select sample IDs
sample_ids <- c("SAMPLE_001", "SAMPLE_002", "SAMPLE_003")

# Multiple sequence alignment
msa_result <- run_msa_from_db(con, sample_ids, method = "Muscle")
alignment_stats(msa_result)

# Build phylogenetic tree
tree <- build_tree_from_db(con, sample_ids, method = "nj")
plot(tree)

# Cluster sequences
clusters <- cluster_from_db(con, sample_ids, threshold = 0.03)
```

### Export Data

```r
# Export to FASTA
export_to_fasta(con, "output.fasta", include_metadata = TRUE)

# Export selected sequences
export_to_fasta(con, "selected.fasta", sample_ids = c("SAMPLE_001", "SAMPLE_002"))

# Export metadata CSV
export_metadata_csv(con, "metadata_export.csv")
```

### Close Connection

```r
close_db_connection(con)
```

## Database Schema

### sequences
| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| sample_id | TEXT | Unique sample identifier |
| sequence | TEXT | DNA sequence |
| sequence_length | INTEGER | Calculated length |
| created_at | DATETIME | Timestamp |

### metadata
| Column | Type | Description |
|--------|------|-------------|
| sample_id | TEXT | Foreign key to sequences |
| sampling_date | DATE | Collection date (YYYY-MM-DD) |
| geo_location | TEXT | Location name |
| geo_county | TEXT | County/region |
| geo_country | TEXT | Country (default: Norway) |
| virus_name | TEXT | Virus name |
| submitter | TEXT | Who submitted |
| sequencing_run | TEXT | Run identifier |
| notes | TEXT | Additional notes |

### analysis_results
| Column | Type | Description |
|--------|------|-------------|
| analysis_name | TEXT | User-defined name |
| analysis_type | TEXT | alignment/phylogeny/clustering |
| sample_ids | TEXT | Comma-separated IDs |
| result_data | TEXT | JSON-serialized results |
| parameters | TEXT | Analysis parameters |

## Moving to Shared Drive

When ready to move to a shared location:

1. Copy the entire `hav_db` folder to the shared drive
2. Update `DB_PATH` in `app.R` (line 19) to point to new location:
   ```r
   DB_PATH <- "N:/shared/hav_db/data/hav.db"
   ```

## Future Development

## Viewer Submodule

The `viewer/` directory is included as a git submodule that points to a fork of the Microreact viewer maintained alongside this project. This lets you keep the viewer's history separate while allowing local edits to be pushed to your fork.

- Fork URL used for this project: `git@github.com:jonbra/viewer.git`
- Topic branch used for local integration changes: `feat/hav-db-integration`

Common commands:

Clone with submodules:
```powershell
git clone --recurse-submodules git@github.com:jonbra/hav_db.git
```

If you've already cloned without submodules:
```powershell
git submodule update --init --recursive
```

Work inside the submodule (make changes, push to your fork):
```powershell
cd viewer
git checkout feat/hav-db-integration
# make edits
git add .
git commit -m "Describe change"
git push origin feat/hav-db-integration

# Back in the superproject, record the new submodule commit
cd ..
git add viewer
git commit -m "Update viewer submodule to latest"
git push
```

To pull updates from upstream (microreact) into your fork:
```powershell
cd viewer
git fetch upstream
git checkout main
git merge upstream/main
git push origin main
```

Notes:
- The superproject stores the exact commit SHA of the submodule. Updating the viewer requires committing the new SHA in the main repo.
- If you prefer not to manage a submodule, the viewer can be vendored into this repo instead (copy files and remove the nested `.git`).


- [ ] AB1 file import using `sangeranalyseR`
- [ ] Quality trimming for raw sequences
- [ ] Reference sequence management
- [ ] BLAST integration
- [ ] Report generation

## License

Internal use - FHI
