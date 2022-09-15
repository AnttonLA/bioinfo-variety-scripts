import os
import sys
import argparse
import pandas as pd
import numpy as np
from liftover import get_lifter

"""
This script creates and fills folders 'bed_files_for_UCSC_browser' and 'significant_hits_only_per_chr' with .bed
files and .tsv tables respectively. It takes a path to a folder where several GWAS output files are stored, and goes
through each of those files looking for hits that are significant (p-value < 5e-8). It the creates, depending on the
request, either .bed files per chr and phenotype, or tables (.txt) with all significant hits per chr.

The idea behind the bed files is that you can easily upload the file to UCSC browser to check all the hits you got
for a hit in a single chromosome. We expect all hits to belong to the same peak or two.

The tables are just summary statistics of the GWAS, separated by chromosome. The idea is that they can be processed
further by dropping duplicates and LDprunning to get a list of unique GWAS hits.

Inputs:
-i: path to the folder containing the GWAS hits (e.g. 'GWAS_hits/')
-p: project prefix (e.g. 'CB_w11_2022')
-o: path to output folder (e.g. 'output/')
--bed: if specified, bed files will be created
--tsv: if specified, tables with significant hits will be created

Outputs:
1. .bed files for UCSC browser (e.g. 'output/bed_files_for_UCSC_browser/')
2. .tsv tables with all significant variants (e.g. 'output/significant_hits_only_per_chr/')
"""

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Path to the folder containing GWAS hits', required=True)
parser.add_argument('-p', '--project', help='Project name. Used to prefix the output file names', required=True)
parser.add_argument('-o', '--output', help='Path to output folder', required=True)
parser.add_argument('--bed', help='Request to make bed files for UCSC browser', action="store_true")
parser.add_argument('--tsv', help='Request to make tsv tables with all significant variants', action="store_true")

args = parser.parse_args()
if not args.bed and not args.tsv:
    print("Please specify at least one output type: --bed or --tsv")
    sys.exit(1)

files_filepath = args.input
project_name = args.project
output_filepath = args.output


def make_bed_file_hg38(df, label: str, out=os.getcwd()):
    """Function that makes bed file out of GWAS results DataFrame, lifting over results to hg38.
    It stores the resulting file in the specified output folder under the name <label>_GWAS_hits_hg38.bed.

    :param df: GWAS results DataFrame
    :param label: label for the bed file
    :param out: path to output folder. Default is current working directory
    :return: None
    """
    filename = out + '/' + label + "_GWAS_hits_hg38.bed"
    converter = get_lifter('hg19', 'hg38')

    with open(filename, 'w') as out_file:
        for index in df.index:

            chromosome_hg37 = 'chr' + str(df.loc[[index]].chromosome.values[0])
            loc_hg37 = int(df.loc[[index]].position.values[0])
            try:
                liftover_tuple_hg38 = converter.convert_coordinate(chromosome_hg37, loc_hg37)[0]
            except IndexError:
                continue
            chrom = liftover_tuple_hg38[0]
            loc = str(liftover_tuple_hg38[1])

            # Get regular rsID from DeCODE rsIDs
            decode_rsid = df.loc[[index]].rsid.values[0]
            if str(decode_rsid).startswith('rs'):
                rsid = str(decode_rsid).split(":")[0]
            else:
                rsid = str(decode_rsid)

            out_file.write(chrom + '\t' + loc + '\t' + '\t' + loc + '\t' + rsid + '\n')  # Add line to bed file


def extract_significant_hits_hg38(in_df, chrom, project_prefix, pval_threshold=5e-8):
    """Function that extracts significant hits from GWAS results DataFrame. It stores the hits into a txt file.

    :param in_df: GWAS results DataFrame
    :param chrom: chromosome number
    :param project_prefix: prefix for the output file name
    :param pval_threshold: p-value threshold. Default is 5e-8
    :return: DataFrame with significant hits
    """
    filename = f'/significant_hits_only_per_chr/{project_prefix}_chr{str(chrom)}_hg38.txt'
    # If file already exists, read it into a df. Otherwise, create empty df.
    if os.path.exists(output_filepath + filename):
        df = pd.read_csv(output_filepath + filename, sep='\t')
        print(f"File {filename} already exists. Appending to it.")
    else:
        df = pd.DataFrame()
        print(f"File {filename} does not exist. Creating it.")

    # Do liftover of in_df to hg38
    converter = get_lifter('hg19', 'hg38')
    in_df.position = in_df.apply(lambda x: converter.convert_coordinate(x.chromosome, x.position)[0][1], axis=1)

    df = pd.concat([df, in_df[in_df.pval < pval_threshold]])  # Concat df read from file and significant hits of in_df
    df.sort_values(by=['chromosome', 'position'], inplace=True)
    # Drop duplicates
    df.drop_duplicates(keep='first', inplace=True)
    df.to_csv(output_filepath + filename, sep='\t', index=False)


# Create folders if they do not exist (and they have been requested)
if not os.path.exists(output_filepath):
    os.makedirs(output_filepath)
if not os.path.exists(output_filepath + '/bed_files_for_UCSC_browser') and args.bed:
    os.mkdir(output_filepath + '/bed_files_for_UCSC_browser')
if not os.path.exists(output_filepath + '/significant_hits_only_per_chr') and args.tsv:
    os.mkdir(output_filepath + '/significant_hits_only_per_chr')

for i, file in enumerate(os.listdir(files_filepath)):
    # Iterate through every file in the specified folder. Ignore the '_mean' MFI files
    if file.endswith('.txt') and not file.endswith('_mean.txt'):
        print(f"Reading ile #{i+1} of {len(os.listdir(files_filepath))}: {file}")
        f = pd.read_csv(files_filepath + '/' + file, sep="\t", chunksize=100000, index_col=False)  # Read by chunks
        list_of_downsampled_chunks = []
        for chunk in f:
            list_of_downsampled_chunks.append(chunk[chunk['pval'] < 5.0e-8].copy())  # P-value threshold
        phenotype = ''.join(file.split('.')[0].split('_')[3:])  # Take phenotype name from file name

        concat_df = pd.concat(list_of_downsampled_chunks)  # Concatenate chunks into a single DataFrame
        concat_df = concat_df.assign(phenotype=phenotype)  # Add phenotype column
        rsid_short_column = pd.Series(rsid.split(':')[0] if rsid.startswith('rs') else rsid for rsid in concat_df.rsid.to_list())
        concat_df.insert(1, 'rsid_short', value=rsid_short_column)
        # concat_df = concat_df[concat_df['pval'] <= 5e-8]  # additional pvalue filter

        chrs_with_hits = list(concat_df['chromosome'].unique())  # list of chrs with p-values that cross threshold
        chrs_with_hits.sort()
        if chrs_with_hits:
            n = concat_df.iloc[0, 9]  # This position of the df should correspond to number of samples of the GWAS
            if type(n) != np.int64:
                raise ValueError(f"Number of samples is not the correct type. The value of n is {n} and its type is {type(n)}")
            for chromosome in chrs_with_hits:
                if args.bed:
                    # Create bed file for UCSC browser
                    make_bed_file_hg38(concat_df[concat_df.chromosome == chromosome].copy(),
                                       label=f'{project_name}_chr{str(chromosome)}_{phenotype}_n{n}_w11_2022',
                                       out=output_filepath + '/bed_files_for_UCSC_browser')
                if args.tsv:
                    # Create table of variants that crossed significance threshold ('hit table')
                    extract_significant_hits_hg38(concat_df[concat_df.chromosome == chromosome].copy(),
                                                  project_prefix=project_name,
                                                  chrom=chromosome)
        else:
            print(f"\tNo hits found")
