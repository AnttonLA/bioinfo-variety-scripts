import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import clipboard


"""
This script creates a plot from the data in the file specified by the user.
The input file should be a .bed-like table with the columns 'MAF' and 'beta', or equivalent.
The script will generate a plot with the Allele Frequencies on the x-axis and the beta values on the y-axis.
"""

parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='The input file to read from')
parser.add_argument('-m', '--maf_column', help='The name of the column containing the MAFs. If none is'
                                               'specified, the default is "MAF".')
parser.add_argument('-b', '--beta_column', help='The name of the column containing the beta values. If none is'
                                                'specified, the default is "beta".')
parser.add_argument('-t', '--title', help='The title of the plot. If none is specified, no title will be added.')
parser.add_argument('-i', '--interactive', action='store_true', help='If this flag is set, the plot will be interactive'
                                                                     ' and the user will be able to click on points to'
                                                                     ' see the corresponding name, MAF and beta values.')
parser.add_argument('-n', '--name_column', help='The name of the column containing the names of the SNPs. If none is'
                                                'specified, the default is "name".')
parser.add_argument('-o', '--output', help='The output file to write to. If none is specified, the plot will be shown '
                                           'on screen, but not saved.')
args = parser.parse_args()

# Check that the input file exists
if not os.path.isfile(args.input_file):
    raise ValueError('The input file does not exist')

# If output name is given but is not a .png file, add the .png extension
if args.output is not None and not args.output.endswith('.png'):
    args.output += '.png'

# Read the input file
df = pd.read_csv(args.input_file, sep='\t')

# Check that the input file contains the required columns
if args.maf_column is None:
    args.maf_column = 'MAF'
if args.beta_column is None:
    args.beta_column = 'beta'
if args.maf_column not in df.columns:
    raise ValueError('The input file does not contain the column "{}"'.format(args.maf_column))
if args.beta_column not in df.columns:
    raise ValueError('The input file does not contain the column "{}"'.format(args.beta_column))

# Add new columns 'MAF_to_plot' and 'beta_to_plot' to the dataframe
df['MAF_to_plot'] = df[args.maf_column]
df['beta_to_plot'] = df[args.beta_column]

# If interactive mode is set, record the names of the SNPs
if args.interactive:
    if args.name_column is None:
        args.name_column = 'name'
    if args.name_column not in df.columns:
        raise ValueError('The input file does not contain the column "{}"'.format(args.name_column))
    df['name'] = df[args.name_column]


# For those entries where MAF > 0.5, flip it, and make the beta negative
df.loc[df[args.maf_column] > 0.5, 'MAF_to_plot'] = 1 - df['MAF_to_plot']
df.loc[df[args.maf_column] > 0.5, 'beta_to_plot'] = -df['beta_to_plot']


sns.set_style('whitegrid')
sns.set_context('paper')
# Create the plot if interactive mode is not set
if not args.interactive:
    fig, ax = plt.subplots()
    ax.scatter(df['MAF_to_plot'], df['beta_to_plot'], s=1)
    ax.set_xlabel('MAF')
    ax.set_ylabel('beta')
    ax.set_xscale('log')
    ax.set_xlim(1e-4, 0.5)

    # Add title if given
    if args.title is not None:
        ax.set_title(args.title)

    # Save or show the plot
    if args.output is None:
        plt.show()
    else:
        fig.savefig(args.output)

# Create the plot if interactive mode is set
else:
    fig = px.scatter(df, x='MAF_to_plot', y='beta_to_plot', hover_data=['name'])

    # define the callback function to copy the clicked point's name to clipboard
    def copy_to_clipboard(trace, points, state):
        name = df.loc[points.point_inds[0], 'name']
        clipboard.copy(name)


    # register the callback function with the interactive plot
    fig.data[0].on_click(copy_to_clipboard)
    # display the interactive plot
    fig.show()
