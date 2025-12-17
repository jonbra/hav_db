Conda environment for HAV DB
---------------------------

This repository can use the `mafft` binary installed from conda for sequence
alignment. To create a dedicated conda environment for this project run:

    conda env create -f conda/environment.yml

Then activate it:

    conda activate hav_db

Confirm `mafft` is available:

    mafft --version

Notes:
- R packages are installed via R (CRAN/Bioconductor) into your R library.
- You do not need to install R packages inside conda; only `mafft` is required
  by the updated alignment code.
