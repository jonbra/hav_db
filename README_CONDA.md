# Conda Environment for HAV DB

This repository uses bioinformatics tools installed via conda for sequence
alignment and phylogenetic analysis. The R code automatically locates the
tools from the `hav_db` conda environment, so you don't need to activate
the environment before running the Shiny app.

## Setup

Create the conda environment:

```bash
conda env create -f conda/environment.yml
```

## Installed Tools

The `hav_db` environment includes:

- **mafft** - Multiple sequence alignment
- **iqtree** (v3.0.1) - Maximum likelihood phylogenetic inference
- **biopython** - Python bioinformatics utilities

## Verify Installation

Activate the environment and check the tools:

```bash
conda activate hav_db
mafft --version
iqtree2 --version
```

## How It Works

The R analysis code (`R/analyses/alignment.R` and `R/analyses/phylogeny.R`)
automatically searches for binaries in the `hav_db` conda environment at
common installation paths:

- `~/miniforge3/envs/hav_db/bin/`
- `~/miniconda3/envs/hav_db/bin/`
- `~/anaconda3/envs/hav_db/bin/`

This means you do **not** need to activate the conda environment before
starting the Shiny app - the tools will be found automatically.

## Updating the Environment

To add new tools or update versions, edit `conda/environment.yml` and run:

```bash
conda env update -f conda/environment.yml --prune
```

## Notes

- R packages are installed via R (CRAN/Bioconductor) into your R library.
- You do not need to install R packages inside conda.
- The conda environment is only used for external bioinformatics tools.
