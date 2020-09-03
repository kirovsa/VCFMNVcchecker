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
Rscript -e "rmarkdown::render('MAFMNVChecker.Rmd',params=list(MAF='testset.short.txt',author='My Name',title='Example format'))"
HLA and immunoglobulin genes are excluded by VCFMNVChecker

Please make sure your MAF file conatins all required columns (consult the included example file).
If columns are missing or mis-named your output file will contain an error message and analysis will fail.

The columns that are mandatory are:

"Symbol","Sample","Transcript","protein_position","HGVSp","t_ref_count","t_alt_count","Codons","Consequence"

The symbol and sample column can be named in a different name as long as they contain the strings "symbol" and "sample". If there are multiple columns containing those strings only the first will be used so make sure it is the relevant one.


