# HAV Database - Copilot Instructions

## Project Overview

HAV Database is a Shiny web application for managing Hepatitis A Virus (HAV) sequences with metadata, supporting phylogenetic analysis and Microreact visualization. It combines an R/Shiny frontend with SQLite storage and external bioinformatics tools (MAFFT, IQ-TREE) from a conda environment.

## Architecture

### Runtime
- **R Shiny app** (`app.R`) on port 3838 - main database/analysis interface
- **Microreact Viewer** (`viewer/`) on port 3000 - *experimental* local viewer for `.microreact` files (minimal functionality; primary workflow uses microreact.org)

### Microreact Integration
The primary visualization workflow uploads projects to **microreact.org** via [R/microreact_api.R](R/microreact_api.R). The local `viewer/` Next.js app is experimental and not actively developed - focus efforts on API integration improvements.

### Key Data Flow
```
FASTA + CSV → parse_metadata.R → cleaned_*.csv → import_database.R → SQLite (data/hav.db)
                                                                            ↓
Shiny UI → db_functions.R ← → analysis_functions.R → MAFFT/IQ-TREE → Microreact API
```

### Core Module Organization
- [R/db_functions.R](R/db_functions.R) - CRUD operations, `get_db_connection()`, `add_sequence_with_metadata()`
- [R/analysis_functions.R](R/analysis_functions.R) - MSA, tree building via external tools
- [R/microreact_api.R](R/microreact_api.R) - API integration for uploading projects to microreact.org
- [R/modules/](R/modules/) - Shiny server modules (analysis_module.R, microreact_module.R)
- [R/analyses/](R/analyses/) - Pure R wrappers: alignment.R (MAFFT), phylogeny.R (NJ/UPGMA/IQ-TREE)

## Development Setup

### Required Steps
```bash
# 1. Create conda environment for bioinformatics tools
conda env create -f conda/environment.yml

# 2. Install R packages (see README.md for full list)
# Key: DBI, RSQLite, shiny, shinydashboard, Biostrings, ape, httr, jsonlite

# 3. Initialize database
Rscript setup_database.R

# 4. Launch both services
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"  # Terminal 1
cd viewer && npm run dev -- -p 3000                    # Terminal 2
```

### Configuration
- Database path: Set via `HAV_DB_PATH` env var or edit `config_local.R` (not committed)
- Default location: `data/hav.db`

## Code Conventions

### External Tool Integration Pattern
R code uses `get_conda_bin()` from [R/analyses/alignment.R](R/analyses/alignment.R) to locate binaries in the `hav_db` conda environment without requiring activation:
```r
bin <- get_conda_bin("mafft")  # Searches ~/miniforge3/envs/hav_db/bin/, etc.
system2(bin, args = c("--auto", input_fasta), stdout = output_fasta)
```

### Database Operations
Always use parameterized queries via DBI:
```r
dbGetQuery(con, "SELECT * FROM sequences WHERE sample_id = ?", params = list(sample_id))
```

### Metadata Import Pipeline
1. Raw exports go to `export/` directory
2. [R/parse_metadata.R](R/parse_metadata.R) normalizes country names (Norwegian→English) and cleans data
3. [import_database.R](import_database.R) loads cleaned CSVs into SQLite

### Shiny Module Pattern
UI tabs in `R/ui_tab_*.R` return `tabItem()` components. Server logic uses `register_*_server()` functions that receive shared `rv` (reactiveValues) and `con` (database connection reactive).

## Database Schema

Four tables in SQLite (see [setup_database.R](setup_database.R)):
- `sequences` - sample_id (unique), sequence, sequence_length
- `metadata` - Foreign key to sequences, includes genotype, geo_country, sampling_date
- `analysis_results` - Stores MSA/tree results as JSON
- `blast_results` - Stores BLASTn hits with SNP counts, links query to hit sequences

## Testing & Validation

- [R/qc_checks.R](R/qc_checks.R) - Database integrity checks: `run_all_checks(db_path)`
- **No formal test suite yet** - use `testthat` framework when adding tests

### Recommended Test Coverage
When building tests, focus on these critical paths:
1. **Database operations** - `add_sequence()`, `add_metadata()`, `get_sequences_by_ids()` round-trip correctly
2. **Import pipelines** - `import_fasta_to_db()` and `import_with_metadata()` handle edge cases (duplicates, malformed FASTA, missing metadata fields)
3. **Analysis functions** - `run_msa()` produces valid alignments, tree builders return valid `phylo` objects
4. **Metadata cleaning** - `parse_metadata.R` country name normalization works for known inputs

Example test structure:
```r
# tests/testthat/test-db_functions.R
test_that("add_sequence stores and retrieves correctly", {

con <- dbConnect(SQLite(), ":memory:")
# setup schema, test, teardown
})
```

## Security Notes

- **Microreact API tokens** start with `eyJ` (JWT format) - never commit to version control
- Store tokens in `config_local.R` (gitignored) or use `HAV_DB_PATH` / `MICROREACT_TOKEN` env vars
- Tokens are session-scoped in the Shiny app; users re-enter on restart

## Key External Dependencies

| Tool | Source | Used For |
|------|--------|----------|
| MAFFT | conda (hav_db env) | Multiple sequence alignment |
| IQ-TREE 3.x | conda (hav_db env) | Maximum likelihood trees |
| BLAST+ | conda (hav_db env) | Sequence similarity search (blastn, makeblastdb) |
| Microreact API | microreact.org | Project upload/sharing |

## Common Tasks

### Adding a New Analysis Type
1. Add pure R function in `R/analyses/`
2. Wire into `R/analysis_functions.R`
3. Add UI controls in `R/ui_tab_analysis.R`
4. Register handlers in `R/modules/analysis_module.R`

### Modifying the Database Schema
1. Update [setup_database.R](setup_database.R)
2. Consider migration strategy (currently drops/recreates tables)
3. Update corresponding functions in `R/db_functions.R`

## Local Deployment Only

This application is designed for local/on-premise use only. There is no remote deployment workflow - run directly on your workstation or local server.
