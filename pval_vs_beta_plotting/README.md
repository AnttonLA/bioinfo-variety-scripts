# pval_vs_beta_plotting

This repository contains scripts for plotting minor allele frequency (MAF) vs beta values for GWAS variants.
The scripts use argparse to parse command line arguments. You can see the arguments by running the script with the -h flag.

## Usage

You can create a conda environment with the required dependencies using the provided environment.yml file:

    conda env create -f environment.yaml

Then activate the environment:

    conda activate pval_vs_beta_plotting

And run the script:

    python create_plotly_iplot.py <input_file> -m <maf_column> -b <beta_column> -n <name_column>\
    -c <chromosome_column> -p <position_column>

## create_plotly_iplot.py

Creates a beta vs MAF plot using plotly. This plot open in your browser and is interactive: if you click on a point it
will open a new tab with the dbSNP page for that SNP and/or* UCSC at that genomic position (it is assumed that your data
is in hg38).

`NOTE that the current version expects MAF values to be percentages, not fractions! (i.e. 50 instead of 0.5)`

\* This has changed several times during development. For now, you can alter the script to suit your needs:)
