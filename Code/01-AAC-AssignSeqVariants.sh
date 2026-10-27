#!/bin/bash
# Script Name: 01-AssignSeqVariants.sh
# Description: This script is a wrapper for mapping single end and paired end data
# Author: E Armstrong
# Date: 2025-04-30
# Useage: sbatch --job-name=AssignSeq --nodes=1 --cpus-per-task=8 --mem-per-cpu=64G --cluster cbsuschmidt /workdir/eja98/BIO5D_Data/01-AAC-AssignSeqVariants.sh
#SBATCH --chdir=/workdir
#SBATCH --output=/workdir/eja98/BIO5D_Data/output_%j.out
#SBATCH --error=/workdir/eja98/BIO5D_Data/error_%j.txt


# Load necessary modules if required by your HPC system
module load R/4.5.2 # Check available R modules on your system

# Run the R script using Rscript
# Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/workdir/eja98/BIO5D_Data/pandoc'); rmarkdown::render('/workdir/eja98/BIO5D_Data/01-AssignSeqVariants.Rmd', 'html_document')"
Rscript /workdir/eja98/BIO5D_Data/01-AAC-AssignSeqVariants.R