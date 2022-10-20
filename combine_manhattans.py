import os
import sys

import pandas as pd
import argparse

"""
This script combines all the GWAS output files in a directory into a single GWAS output file. The output will be
stored in the same directory where the input files, in a folder named 'combined_output/'.
There are three separate steps for the processing:
    1. Extract all the significant hits (below specified p-value threshold) from the GWAS output files.
    2. Combine all the significant hits into a single file.
    3. Take a GWAS output file and replace the p-values with the p-values from the combined file in the relevant rows.
"""

parser = argparse.ArgumentParser(description="Combine GWAS sumstats into one file for a single Manhattan plot.")
parser.add_argument('GWAS_dir', metavar="path",
                    help="Directory where the summary statistics files we want to merge are located", type=str)
parser.add_argument('--prefix', help="Prefix for the output files. Used to differentiate between projects."
                    "For example 'BV' or 'CB'", type=str, default="")
parser.add_argument("-s", "--steps", nargs="+", type=str, help="Numbers of the steps you want to run. Space delimited.")
parser.add_argument("-p", "--pval", type=int,  # TODO : make this into a float instead
                    help="Exponent of the p-value threshold to extract variants from GWAS output files")
parser.add_argument("-r", "--include_repeats", type=bool, default=True,
                    help="whether to include all the associations for a variant, or just the most significant one")
parser.add_argument("-a", "--alias_file", type=str, default=None,
                    help="Path to the file with aliases for the phenotypes (for example, to simplify to lineages)."
                         "If not specified, the script will add a column identical to the 'phenotype' column to the"
                         "final table.")
# Immune GWAS alias file: '/home/antton/Projects/Immune_GWAS/data/processed/phenotype_lineage_curated.csv'
# TODO : make sure this is correct

args = parser.parse_args()


files_filepath = args.GWAS_dir  # The directory containing the GWAS output files to combine.
project_id = args.prefix  # Prefix used for file names. For example "BV" or "CB"
steps_to_run = args.steps  # The steps to run. For example "1" or "123"

print("\nGiven filepath: ", files_filepath)
print("Project prefix: ", project_id)
print("Steps requested: ", steps_to_run)

# Make sure 'steps_to_run' is a string and contains only numbers from 1 to 3
for step in steps_to_run:
    if step not in ['1', '2', '3']:
        print("ERROR: Step numbers must be 1, 2 or 3. You entered: ", steps_to_run)
        sys.exit(1)

p_value_threshold = args.pval  # Exponent of the p-value threshold. For example "5" for 1e-5.

# Make sure the inputted p-value threshold is a valid number.  # TODO : make this into a float instead
try:
    int(p_value_threshold)  # TODO: This typecheck is no longer needed since argparse takes care of it
except ValueError:
    print("The p-value threshold must be an integer.")
    sys.exit(1)

p_value_threshold = 10**-(int(p_value_threshold))  # Actual p-value threshold, not the exponent.

print("P-value threshold: ", p_value_threshold)

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
    print('\nCommencing STEP 1: filling all_significant_variant_tables folder...')
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
    print('\nCommencing STEP 2: Combining all tsv files...')

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

    if args.include_repeats:  # each variant only once
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
    print('\nCommencing STEP 3: Replacing template GWAS output file with lowest p-values...')

    hits_file = output_filepath + '/' + hit_table_file_name  # TODO: If I use this variable then STEP 2 must be run :(
    template_manhattan_file = output_filepath + '/template_manhattan.txt'
    output_manhattan_file = output_filepath + '/output_manhattan.txt'

    # Check if no template GWAS file exist
    if not os.path.isfile(template_manhattan_file):
        # Create template GWAS file by copying the last file in os.listdir(files_filepath) to template_manhattan_file
        print("No template GWAS file found. Creating one...")
        with open(template_manhattan_file, 'w') as f:
            if os.listdir(files_filepath)[-1].endswith('.txt'):
                f.write(open(files_filepath + '/' + os.listdir(files_filepath)[-1]).read())
            else:
                # Throw error if os.listdir(files_filepath)[-1] is not a .txt file
                raise ValueError('The last file in the directory is not a .txt file.')
                # TODO: looking for the last file in the directory is completely arbitrary (and hacky). It relies on the
                #  file structure being kept. A better way must exist to find a GWAS file.
        print("File created.  Replacing rows...")

    else:
        print("Template GWAS file found. Replacing rows...")

    hits_df = pd.read_csv(hits_file, sep='\t', index_col=None)
    hits_df.drop(columns=['rsid_short'], inplace=True)  # Make it identical to the GWAS_df

    # If there is an alias file, read it and add the alias column to the hits_df
    if args.alias_file:
        if not os.path.isfile(args.alias_file):
            raise ValueError('The provided alias file does not exist.')

        # Dictionary used to assign alias to each phenotype
        alias_dict = {}
        with open(args.alias_file, 'r') as in_file:
            for line in in_file:
                split_line = line.split(',')
                alias_dict[split_line[1].replace(' ', '')] = split_line[2].rstrip()

        # add column 'alias' to the hits_df dataframe using the alias_dict
        hits_df['alias'] = hits_df['phenotype'].map(alias_dict)

    else:
        hits_df['alias'] = hits_df['phenotype']

    # Read in the template GWAS file and replace the relevant variant entries with the contents of the hits_df dataframe
    f = pd.read_csv(template_manhattan_file, sep="\t", chunksize=100000)
    list_of_downsampled_chunks = []
    for chunk in f:
        list_of_downsampled_chunks.append(chunk.copy())
    templateGWAS_df = pd.concat(list_of_downsampled_chunks)

    # Replace the relevant entries in the templateGWAS_df dataframe with the contents of the hits_df dataframe
    out_df = pd.merge(templateGWAS_df, hits_df, how='outer',
                      on=["rsid", "chromosome", "position", "A1", "A2", "pval", "beta", "tstat", "n"])
    # TODO: This merge is not working properly. It is not replacing the rows in the templateGWAS_df dataframe
    # Write the output_manhattan_file
    out_df.to_csv(output_manhattan_file, sep='\t', index=False)

    # You can check the rows where the phenotype and alias are not NaN with the bash command:
    # awk '$11' output_manhattan.txt

    print("STEP 3 done! Output file: ", output_manhattan_file)
