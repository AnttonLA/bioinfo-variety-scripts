import os
import sys

import pandas as pd
import numpy as np

"""
This script combines all the GWAS output files in a directory into a single GWAS output file. The output will be
stored in the same directory where the input files, in a folder named 'combined_output/'.
There are three separate steps for the processing:
    1. Extract all the significant hits (below specified p-value threshold) from the GWAS output files.
    2. Combine all the significant hits into a single file.
    3. Take a GWAS output file and replace the p-values with the p-values from the combined file in the relevant rows.
"""

args = sys.argv[1:]
if len(args) != 5:
    print("Usage: combine_manhattans.py <GWAS data directory> <project_prefix> <steps to run> <p-value threshold>"
          " <drop_repeats boolean>")
    sys.exit(1)

files_filepath = args[0]  # The directory containing the GWAS output files to combine.
project_id = args[1]  # Prefix used for file names. For example "BV" or "CB"

steps_to_run = args[2]  # The steps to run. For example "1" or "123"

# Make sure 'steps_to_run' is a string and contains only numbers from 1 to 3
if not steps_to_run.isdigit() or not (1 <= int(steps_to_run[0]) <= 3):  # Not the best handling, but will do for now.
    print("'steps_to_run' must be a string and contain only numbers from 1 to 3")
    sys.exit(1)

p_value_threshold = args[3]  # Exponent of the p-value threshold. For example "5" for 1e-5.

# Make sure the inputted p-value threshold is a valid number.
try:
    int(p_value_threshold)
except ValueError:
    print("The p-value threshold must be an integer.")
    sys.exit(1)

p_value_threshold = 10**-(int(p_value_threshold))  # Actual p-value threshold, not the exponent.

# Whether to drop repeated hits. For example "True" or "False". Default is "True".
if args[4].lower() in ["true", "t", "yes", "y"]:
    drop_repeats = True
elif args[4].lower() in ["false", "f", "no", "n"]:
    drop_repeats = False
else:
    print("The drop_repeats argument must be either 'True' or 'False'.")
    sys.exit(1)

if files_filepath.endswith("/"):
    files_filepath = files_filepath[:-1]
output_filepath = files_filepath + '/combined_output'

# Create 'combined_output' folder if it doesn't exist
if not os.path.exists(output_filepath):
    os.makedirs(output_filepath)
# Create 'all_significant_variant_tables' folder if it doesn't exist
if not os.path.exists(output_filepath + '/all_significant_variant_tables'):
    os.mkdir(output_filepath + '/all_significant_variant_tables')

########################################################################################################################
# STEP 1
########################################################################################################################
# This step fills the folder 'all_significant_variant_tables' with .tsv tables.
# Out of each GWAS output file, it only takes the entries that have a p-value above the inputted threshold.
# For each phenotype, each chromosome is outputted into a separate file.

# Inputs: GWAS output files

# Outputs: .tsv files in the folder 'all_significant_variant_tables'
week_year_str = "w11_2022"

if '1' in steps_to_run:
    print('Commencing STEP 1: filling all_significant_variant_tables folder...')
    for i, gwas_file in enumerate(os.listdir(files_filepath)):
        if gwas_file.endswith('.txt'):
            print(f"File {i + 1} of {len(os.listdir(files_filepath))}")
            f = pd.read_csv(files_filepath + '/' + gwas_file, sep="\t", chunksize=100000)
            list_of_downsampled_chunks = []
            for chunk in f:
                list_of_downsampled_chunks.append(chunk[chunk['pval'] < p_value_threshold].copy())
                # Note the p-value threshold
            phenotype = ''.join(gwas_file.rstrip('.txt').split('_')[3:])
            concat_df = pd.concat(list_of_downsampled_chunks)
            concat_df = concat_df.assign(phenotype=phenotype)
            # Add shortened rsID
            concat_df.insert(1, 'rsid_short',
                             [rsid.split(':')[0] if rsid.startswith('rs') else rsid for rsid in concat_df.rsid.to_list()])
            # concat_df = concat_df[concat_df['pval'] <= 5e-8]  # additional pvalue filter if needed
            chrs_with_hits = list(
                concat_df['chromosome'].unique())  # list of chromosomes with p-values that cross threshold
            chrs_with_hits.sort()
            if chrs_with_hits:
                n = concat_df.iloc[0, 9]  # get N from dataframe
                for chromosome in chrs_with_hits:
                    concat_df[concat_df.chromosome == chromosome].to_csv(output_filepath
                                                                         + f'/all_significant_variant_tables/{project_id}_GWAS_chr{str(chromosome)}_{phenotype}_n{n}_{week_year_str}_hit_table.tsv',
                                                                         sep='\t', index=False)

########################################################################################################################
# STEP 2
########################################################################################################################
# This step goes over every .tsv file in the 'all_significant_variant_tables' folder and concatenates all into a df.
# This df is then sorted by p-value, and for each distinct rsID only the lowest p-value is kept.

# Inputs: .tsv files in the folder 'all_significant_variant_tables/'

# Outputs: "hit_table_file" file (name depends on the 'drop_repeats' boolean)

if '2' in steps_to_run:
    print('Commencing STEP 2: Combining all tsv files...')

    full_hits_table_df = pd.DataFrame()

    # Check if directory 'all_significant_variant_tables' exists
    if not os.path.isdir(output_filepath + '/all_significant_variant_tables'):
        print('ERROR: all_significant_variant_tables folder does not exist.')
        print('Please run step 1 first so that the directory is created and populated.')
        sys.exit()

    path_to_tsv_files_folder = output_filepath + '/all_significant_variant_tables'

    # Iterate through the files, concat contents into a df
    num_all_tables = len(os.listdir(path_to_tsv_files_folder))
    for i, hit_table_file in enumerate(os.listdir(path_to_tsv_files_folder)):
        if hit_table_file.endswith('.tsv'):

            # Table made from current file
            table_df = pd.read_csv(path_to_tsv_files_folder + '/' + hit_table_file, sep="\t", index_col=0)
            full_hits_table_df = pd.concat([full_hits_table_df, table_df], sort=False)

            if num_all_tables > 100 & i % 100 == 0:  # Print progress. Every 100 files >= 100 files, otherwise every 10)
                print(f"Table  {i} / {num_all_tables}")
            elif num_all_tables > 10 & i % 10 == 0:
                print(f"Table  {i} / {num_all_tables}")

            # if i >300:  # for testing
            #   break

    if drop_repeats:  # each variant only once
        print("Dropping repeated rsid-s. Keeping only most significant association per variant")
        # for each row with the same rsid, keep the one with the lowest p-value
        full_hits_table_df = full_hits_table_df.sort_values(by=['pval']).groupby('rsid').head(1)
        # Note that since we only take one row, we fail to gather all the different phenotypes where the variant is
        # significant. This is not necessarily true, since a variant can have associations with several traits.
        # In this version of the final table, each row belongs to a distinct rsID, showing the stats of the phenotype
        # that gave the lowest p-value for that SNP.

        full_hits_table_df = full_hits_table_df.sort_values(by=['chromosome', 'position'])  # sort by column chromosome

        hit_table_file_name = f'{project_id}_GWAS_hits_only_10E{p_value_threshold}_{week_year_str}'\
                              '_one_pheno_only_per_variant.txt '
        print("Rows of the final table containing all our hits (counting each variant only once): ",
              len(full_hits_table_df.index.to_list()))

    else:  # each variant can appear multiple times if it has multiple associations with different phenotypes
        print("Keeping all association for each variant")
        full_hits_table_df = full_hits_table_df.sort_values(by=['chromosome', 'position'])  # sort by column chromosome
        hit_table_file_name = f'{project_id}_GWAS_hits_only_10E{p_value_threshold}_{week_year_str}'\
                              '_all_phenotypes_per_variant.txt '
        print("Rows of the final table containing all our hits (same variant can be counted several times): ",
              len(full_hits_table_df.index.to_list()))

    # Save the final table to file
    full_hits_table_df.to_csv(output_filepath + '/' + hit_table_file_name, sep="\t")


########################################################################################################################
# STEP 3
########################################################################################################################
# This step takes a GWAS output file and replaces some rows with the rows from the file named in the
# 'hit_table_file_name' variable

# Inputs: template GWAS output file, file named in variant 'hit_table_file_name'

# Outputs: GWAS output file with the rows from the "hit_table_file_name" file replaced into it

if '3' in steps_to_run:
    print('Commencing STEP 3: Replacing template GWAS output file with lowest p-values...')

    hits_file = output_filepath + '/' + hit_table_file_name
    template_manhattan_file = output_filepath + '/template_manhattan.txt'
    output_manhattan_file = output_filepath + '/output_manhattan.txt'

    # Check if no template GWAS file exist
    if not os.path.isfile(template_manhattan_file):
        # Create template GWAS file by copying the last file in os.listdir(files_filepath) to template_manhattan_file
        with open(template_manhattan_file, 'w') as f:
            if os.listdir(files_filepath)[-1].endswith('.txt'):
                f.write(open(files_filepath + '/' + os.listdir(files_filepath)[-1]).read())
            else:
                # Throw error if os.listdir(files_filepath)[-1] is not a .txt file
                raise ValueError('The last file in the directory is not a .txt file.')
                # TODO: looking for the last file in the directory is completely arbitrary (and hacky). It relies on the
                #  file structure being kept. A better way must exist to find a GWAS file.

    # Dictionary used to assign lineage group to each phenotype
    lineage_dict = {}
    with open('/home/antton/Projects/Immune_GWAS/data/processed/phenotype_lineage_curated.csv', 'r') as in_file:
        for line in in_file:
            split_line = line.split(',')
            lineage_dict[split_line[1].replace(' ', '')] = split_line[2].rstrip()

    hits_df = pd.read_csv(hits_file, sep='\t')
    hits_df.drop(columns=['rsid_short'], inplace=True)  # Make it identical to the GWAS_df
    hits_df.rename(columns={'Phenotype': 'phenotype'}, inplace=True)  # Make column name lowercase
    # add column 'lineage' to the hits_df dataframe using the lineage_dict
    hits_df['lineage'] = hits_df['phenotype'].map(lineage_dict)

    # Read in the template GWAS file and replace the relevant variant entries with the contents of the hits_df dataframe
    f = pd.read_csv(template_manhattan_file, sep="\t", chunksize=100000)
    list_of_downsampled_chunks = []
    for chunk in f:
        list_of_downsampled_chunks.append(chunk.copy())
    templateGWAS_df = pd.concat(list_of_downsampled_chunks)

    # Replace the relevant entries in the templateGWAS_df dataframe with the contents of the hits_df dataframe
    out_df = pd.merge(templateGWAS_df, hits_df, how='outer',
                      on=["rsid", "chromosome", "position", "A1", "A2", "pval", "beta", "tstat", "n"])

    # Write the output_manhattan_file
    out_df.to_csv(output_manhattan_file, sep='\t', index=False)

    # You can check the rows where the phenotype and lineage are not NaN with the bash command:
    # awk '$11' output_manhattan.txt

    print("STEP 3 done! Output file: ", output_manhattan_file)
