# Pipeline to cluster variants and assign closest gene

The following is a pipeline to create a "hits table" out of GWAS results. This is a table that describes the suspected **lead variant**, and it's putative **target gene**, also giving summary GWAS statistics for the variant.

For any number of GWASs, it takes the output GWAS tables, takes only the variants that cross the significance threshold, and separates them into chromosomes.
It then assigns a closest gene to each GWAS hit, and compiles the hits into 'peaks' that should correspond to the peaks in a Manhattan plot.
It then produces yet another table by compiling all of the 'peaks' that share common closest genes.

### 1-Get .bed files and full tables for all SNPs above significant threshold, separated by phenotype and chromosome
The script ``generate_bed_files_and_tables.py`` accomplishes this.

### 2-Assign the closest gene to each of the variants in each of the files. 
Also, cluster the variants by distance. Variants withing 100kbs of one another get clustered.
The script ``closest_gene.sh`` does this.
It creates two folders: ``./bed_files_with_closest_gene/`` and ``./bed_files_clustered/``

### 3-For each cluster, take the variant with the lowest p-value ("lead variant")
List all the potentially associated genes.
Produce a final table with the lead variants only, each showing all the associated phenotypes and potential causal genes.
The script ``generate_all_hits_table.py`` generates a file ``all_hits_table.tsv`` containing every SNP from the clustered
.bed files.
The script ``generate_unique_hits_table.py`` then takes this table, and collapses the entries pertaining to the same SNPs
Finally, the script ``get_actual_unique_hits_group_by_genes.py`` will collapse this table even further, by allowing only
one entry for each gene/group of genes associated to the variants.


The script ``run_whole_pipeline.sh`` runs all of these steps one after the other. This will be replaced by a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow in the future.

---

# Individual scripts:

## generate_bed_files_and_tables.py
**Description:** Create and fill folders ``bed_files_for_UCSC_browser/`` and ``all_significant_variant_tables/`` with .bed files and .tsv tables respectively. 

**Language:** Python

**Usage:** ``generate_bed_files_and_tables.py <path to GWAS hits> <project_name> <path to output folder>``

## closest_gene.sh
**Description:** Take a collection of .bed files and find the closest gene for each location. Also cluster the entries by proximity. Requires [samtools](http://www.htslib.org/).

**Language:** Bash

**Usage:** ./closest_gene.sh

## generate_all_hits_table.py
**Description:** Generate a file ``all_hits_table.tsv`` containing every SNP from the clustered .bed files.

**Language:** Python

**Usage:** ``generate_all_hits_table.py <clustered_bed_dir> <closest_gene_dir> <hit_table_dir> <output_dir>``

## generate_unique_hits_table.py
**Description:**  Take ``all_hits_table.tsv`` and compile every entry pertaining to the same variant into a single row, storing all the phenotypes where it crossed significance, as well as all the possible target genes across all GWASs. Put the output in a file named ``unique_hits_table.tsv``.

**Language:** Python

**Usage:** ``generate_unique_hits_table.py <all hits table file> <output_dir>``

## get_actual_unique_hits_group_by_genes.py
**Description:** Take the table ``unique_hits_table.tsv`` and compile it even further by taking all the rows that at least partially share putative target genes. Then take the entry with lowest p-value out of these an save it into a table. The script also takes a p-value cutoff number, so that only the entries that cross that threshold will be considered.

**Language:** Python

**Usage:** ``generate_unique_hits_table.py <path_to_hits_table> <pval_threshold> <output_file>``

## run_whole_pipeline.sh
**Description:** Run the entire hit table generating pipeline at once.

**Language:** Bash

**Usage:** ``./run_whole_pipeline.sh``

