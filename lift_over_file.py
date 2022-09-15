import sys
import os
import pandas as pd
from liftover import get_lifter
import numpy as np

"""LiftOver the input file from hg19 to hg38.
There are probably FAR more efficient and elegant ways to do this, like using the UCSC LiftOver tool
(https://genome.ucsc.edu/cgi-bin/hgLiftOver). This is just a quick and dirty way to move hg19 positions to hg38.
It uses pandas and the liftover package and it is rather slow :)
"""

if len(sys.argv) != 5:
    raise ValueError('Wrong input. Expected usage: lift_over_file.py input_file_path chr_col_name pos_col_name '
                     'output_file_path')

path_to_file = sys.argv[1]
chr_col_name = sys.argv[2]
pos_col_name = sys.argv[3]
output_filename = sys.argv[4]

# Make sure the file exists
if not os.path.isfile(path_to_file):
    raise ValueError('File does not exist: ' + path_to_file)

df = pd.read_csv(path_to_file, sep='\t')
original_column_order_list = df.columns.tolist()  # Save original column order for later

print(f"File loaded. {len(df.index)} rows found. Starting lift over.")
# Make sure chr_col_name and pos_col_name are in the dataframe
if chr_col_name not in df.columns:
    raise ValueError('Column does not exist: ' + chr_col_name)
if pos_col_name not in df.columns:
    raise ValueError('Column does not exist: ' + pos_col_name)

converter = get_lifter('hg19', 'hg38')


def lift_over(chrom, pos):
    """Get liftover for a single position. Returns position only, no chromosome!
    """
    converter_out = converter.query(chrom, pos)
    return converter_out


# This is a bit tricky, it takes a few steps:
# Liftover and get a list with a tuple inside with the new position
df['hg38_pos_full_LO'] = df.apply(lambda x: lift_over(x[chr_col_name], x[pos_col_name]), axis=1)

# Extract only the position
df['hg38_pos'] = df['hg38_pos_full_LO'].str[0].str[1]

# Cast to int after converting NaNs to 0
df['hg38_pos'] = df['hg38_pos'].fillna(0).astype(int)
#df['hg38_pos'].replace(0, np.nan)  # Replace 0s back to Nan
df = df[df['hg38_pos'] != 0]  # Remove rows with 0 as position

# Drop old hg37 position and rename new one to just position. Reorder columns to og order.
df.drop(pos_col_name, axis=1, inplace=True)
df.rename({'hg38_pos': pos_col_name}, axis=1, inplace=True)
df = df.reindex(columns=original_column_order_list)

df.sort_values([chr_col_name, pos_col_name], inplace=True)  # Sort by chr, pos. TODO: Make this optional
print("Lftover complete. Writing output to: ", output_filename)
df.to_csv(output_filename, sep='\t', index=None)
