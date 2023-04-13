# pval_vs_beta_plotting

This repository contains scripts for plotting minor allele frequency (MAF) vs beta values for GWAS variants.
The scripts use argparse to parse command line arguments. You can see the arguments by running the script with the -h flag.

## create_seaborn_plot.py

Creates a beta vs MAF plot using seaborn. The "interactive" flag -i plots the same plot in your browser using plotly,
but no actual interaction has been implemented. Instead, use the create_plotly_iplot.py script for an interactive plot.

## create_plotly_iplot.py

Creates a beta vs MAF plot using plotly. This plot open in your browser and is interactive: if you click on a point it
will open a new tab with the dbSNP page for that SNP.