import os
import polars as pl
import argparse
from liftover import get_lifter

"""
LiftOver the input file from hg19 to hg38.
Specify the column names for chromosome and position (in hg19) amd it will be converted to hg38.
Uses the liftover package (https://github.com/jeremymcrae/liftover) and polars (https://github.com/pola-rs/polars).

You might also want to check out the online UCSC LiftOver tool (https://genome.ucsc.edu/cgi-bin/hgLiftOver).
"""

parser = argparse.ArgumentParser(description='LiftOver the input file from hg19 to hg38. Provide a tab separated file'
                                             'with a column for chromosome and a column for position in hg19.')
parser.add_argument('input_file_path', type=str, help='Path to input file')
parser.add_argument('-c', '--chr_col_name', type=str, required=True, help='Name of the chromosome column')
parser.add_argument('-p', '--pos_col_name', type=str, required=True, help='Name of the position column')
parser.add_argument('-n', '--new_pos_col_name', type=str, help='Name to replace the position column with. '
                                                               'OPTIONAL. If not specified, the output will keep the '
                                                               'name of the original position column.')
parser.add_argument('-s', '--sort', action='store_true', help='Sort output by chromosome and position. '
                                                              'OPTIONAL. If this flag is not used, the output will '
                                                              'keep the original row order.')
parser.add_argument('-k', '--keep_failed', action='store_true', help='Keep rows that failed to liftover. '
                                                                     'OPTIONAL. If this flag is not used, rows that '
                                                                     'failed to liftover will be dropped.')
parser.add_argument('-o', '--output_file_path', type=str, required=True, help='Path to output file')
args = parser.parse_args()

# Check if input file exists
if not os.path.isfile(args.input_file_path):
    raise ValueError('File does not exist: ' + args.input_file_path)

df = pl.read_csv(args.input_file_path, separator='\t', has_header=True)
original_column_order_list = df.columns # Save original column order for later

print(f"File loaded. {df.shape[0]} rows found. Starting lift over.")

# Make sure args.chr_col_name and args.pos_col_name are in the dataframe
if args.chr_col_name not in df.columns:
    raise ValueError('Column does not exist: ' + args.chr_col_name)
if args.pos_col_name not in df.columns:
    raise ValueError('Column does not exist: ' + args.pos_col_name)

# TODO: check column type. It should be able to handle chromosome columns with chr prefix.

converter = get_lifter('hg19', 'hg38')

chr_pos_df = df.select([args.chr_col_name, args.pos_col_name])


def converter_wrapper(chrom, pos):
    """
    Wrapper for the liftover converter. This function is used to apply to the rows of the dataframe.
    Necessary to circumvent failed liftover entries and return None in those cases.

    :param chrom:
    :param pos:
    :return:
    """
    liftover_output = converter[chrom][pos]
    if liftover_output is not None and len(liftover_output) > 0:
        return liftover_output[0][1]
    else:
        return None


# Do the liftover. Add new column with hg38 position with a temporary name 'pos_hg38'
df = df.with_columns(chr_pos_df.map_rows(lambda t: converter_wrapper(t[0], t[1])).to_series().alias("pos_hg38"))

# Print some statistics. Percentage of entries of pos_hg38 that are None
null_count = df.get_column("pos_hg38").null_count()
all_entries_count = df.shape[0]

if null_count > 0:
    print(f"{null_count} entries out of {all_entries_count} could not be lifted over. "
          f"({str(round(null_count/all_entries_count*100, 2))}%)")

    if not args.keep_failed:
        # Drop rows that failed to liftover
        df = df.drop_nulls(subset=["pos_hg38"])
        print(f"Dropped {all_entries_count - df.shape[0]} rows that failed to liftover.")
    else:
        print(f"Rows that failed to liftover were kept with null position value.")

else:
    print("All entries were lifted over successfully.")

# Drop old hg19 position column.
df = df.drop(args.pos_col_name)

# Rename hg38 position column to either the original name or the user specified name.
if args.new_pos_col_name is not None:
    df = df.rename({'pos_hg38': args.new_pos_col_name})
    # Replace the name of the old position column with the new name in the original_column_order_list
    original_column_order_list[original_column_order_list.index(args.pos_col_name)] = args.new_pos_col_name
else:
    df = df.rename({'pos_hg38': args.pos_col_name})

# Make sure the output has the same column order as the input had
df = df.select(original_column_order_list)

if args.sort:
    df = df.sort([args.chr_col_name, args.pos_col_name])  # Sort by chr, pos.

print("Liftover complete!\nWriting output to: ", args.output_file_path)
df.write_csv(args.output_file_path, separator='\t', include_header=True)
