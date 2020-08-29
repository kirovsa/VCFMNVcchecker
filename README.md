# VCFMNVcchecker

Rshiny app intended to check for potential MNV annotation errors and R knitr application that can merge suspected SNVs into MNVs

Requirements:
For the Rshiny tool (vcfMNVCheckRS.R) please install ggplot2, dplyr, stringr, shinyjs and shinyWidgest
For knitr package install ggplot2, Biostrings and kableExtra

vcfMNVCheckRS expects a regular, annotated VCF file with the following fields present:
DP and AD or AF or FA
in the case of strelka, AU, GU, CU, TU should all be present
In the annotation the actual protein substitution is expected, such as p.Ala109Gly

For MAF like files, follow the enclosed file example from TCGA (downloaded from CBio)

To run VCFMNVChecker from command line:
Rscript -e "rmarkdown::render('$PATHTOSCRIPTR/VCCFMNVChecker.Rmd',params=list(MAF='$PATHTODATA/tcgatest.txt`))"

HLA and immunoglobulin genes are excluded by VCFMNVChecker
