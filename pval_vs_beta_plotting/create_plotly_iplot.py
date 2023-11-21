import os
import sys
import argparse
import pandas as pd
import dash
from dash import dcc, html
import plotly.graph_objs as go
import webbrowser

"""
This script creates an interactive plotly scatter plot from a BED-like tab-separated file containing information about
the MAFs, beta values and position of SNPs. It can also take variant name and p-value information.
The script visualizes the variants in a scatter plot, where the x-axis represents the MAFs and the y-axis represents the
beta values. The user can hover over any of the points to see more information about the variant. Clicking on a point
will open the dbSNP page corresponding to the variant. 
"""

parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='The input file to read from')
parser.add_argument('-m', '--maf_column',
                    help='The name of the column containing the MAFs. If none is specified, the default is "MAF".',
                    default='MAF')
parser.add_argument('-b', '--beta_column',
                    help='The name of the column containing the beta values. If none is specified, the default is '
                         '"beta".',
                    default='beta')
parser.add_argument('-n', '--name_column',
                    help='The name of the column containing the names of the SNPs. If none is specified, the default '
                         'is "name".',
                    default='name')
parser.add_argument('-c', '--chromosome_column',
                    help='The name of the column containing the chromosome in which the SNP is located. If none is '
                         'specified, the default name is "chromosome".',
                    default='chromosome')
parser.add_argument('-p', '--position_column',
                    help='The name of the column containing the position of the SNP in the chromosome. If none is '
                         'specified, the default name is "position".',
                    default='position')
parser.add_argument('-pval', '--p_value_column',
                    help='The name of the column containing the p-values of the SNPs. Only used if specified.',
                    default=None)
parser.add_argument('-gene', '--gene_column',
                    help='The name of the column containing the gene names. Only used if specified.',
                    default=None)
args = parser.parse_args()

# Check that the input file exists
if not os.path.isfile(args.input_file):
    raise ValueError('The input file does not exist')

# Read the input file
df = pd.read_csv(args.input_file, sep='\t')


def check_df_columns():
    """
    Check that the input file contains the columns specified by the user or the default column names.
    :return:
    """
    for argument in [args.maf_column, args.beta_column, args.chromosome_column, args.position_column]:
        if argument not in df.columns:
            raise ValueError('The input file does not contain the column "{}"'.format(argument))
    if args.p_value_column and args.p_value_column not in df.columns:
        raise ValueError('The input file does not contain the column "{}"'.format(args.p_value_column))
    if args.gene_column and args.gene_column not in df.columns:
        raise ValueError('The input file does not contain the column "{}"'.format(args.gene_column))


check_df_columns()

# If the entries in 'chromosome' start with 'chr', remove it
if type(df[args.chromosome_column][1]) == str and df[args.chromosome_column][1].startswith('chr'):
    df[args.chromosome_column] = df[args.chromosome_column].str.replace('chr', '')

# If no name column was given, create one by concatenating the chromosome and position
if args.name_column not in df.columns:
    df[args.name_column] = df[args.chromosome_column].astype(str) + ':' + df[args.position_column].astype(str)
    sys.stdout.write('No name column with name "{}" could be found. Variants will be given a name based on the '
                     'chromosome and position columns instead.\n'.format(args.name_column))
else:
    df[args.name_column] = df[args.name_column].astype(str) + ' (' + df[args.chromosome_column].astype(str) + ':' + \
                           df[args.position_column].astype(str) + ')'  # Add the chromosome and position to the name

# Add new columns to the dataframe that will be used for plotting. This way we also keep the original data untouched.
df['MAF_to_plot'] = df[args.maf_column].copy()
df['beta_to_plot'] = df[args.beta_column].copy()
df['name'] = df[args.name_column].copy()
df['chromosome'] = df[args.chromosome_column].copy()
df['position'] = df[args.position_column].copy()
if args.p_value_column:
    df['p_value'] = df[args.p_value_column].copy()
else:
    df['p_value'] = None
if args.gene_column:
    df['gene'] = df[args.gene_column].copy()
else:
    df['gene'] = None


# For those entries where MAF > 0.5, flip it, and make change the sign of the beta
# TODO: be mindful that data can come both in percentages and fractions!!! (e.g. 0.5 vs 50). This should be harmonized!
mask = df[args.maf_column] > 50
df.loc[mask, 'beta_to_plot'] = -df.loc[mask, 'beta_to_plot']
df.loc[mask, 'MAF_to_plot'] = 100 - df.loc[mask, 'MAF_to_plot']


def open_dbsnp_page(chromosome, position):
    """
    Open the dbSNP page for the given chromosome and position
    """
    url = f"https://www.ncbi.nlm.nih.gov/snp/?term={chromosome}%5BChromosome%5D+AND+{position}%5BCHRPOS%5D"
    webbrowser.open_new_tab(url)


def open_ucsc_page(chromosome, position):
    """
    Open the UCSC page for the given chromosome and position
    :param chromosome:
    :param position:
    :return:
    """
    zoom_out = 10_000
    start_zoom = max(0, int(position) - zoom_out)
    end_zoom = int(position) + zoom_out
#    url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr{chromosome}%3A{position}-{position}"
    url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr{chromosome}%3A{start_zoom}-{end_zoom}&highlight=chr{chromosome}%3A{position}-{position}"
    webbrowser.open_new_tab(url)


# Create the Dash app
app = dash.Dash(__name__)

# Create the scatter plot
scatter = go.Scatter(
    x=df['MAF_to_plot'],
    y=df['beta_to_plot'],
    mode='markers',
    marker=dict(
        size=15,
        color='blue',
        opacity=0.5,
    ),
    text=df['name'],
    customdata=df['p_value'],
)

fig = go.Figure(scatter)
if args.p_value_column:
    hovertemplate_string = 'variant: %{text}<br>MAF: %{x}<br>β: %{y}<br>p-value: %{customdata}'
else:
    hovertemplate_string = 'variant: %{text}<br>MAF: %{x}<br>β: %{y}'
fig.update_traces(hovertemplate=hovertemplate_string)


# Define the callback function
@app.callback(
    dash.dependencies.Output('scatter', 'figure'),
    [dash.dependencies.Input('scatter', 'clickData')],
    [dash.dependencies.State('scatter', 'figure')]
)
def update_point(clickData, figure):
    """
    Update the point that was clicked on, making it red and small. Also open the dbSNP page for the variant.
    :param clickData: Dictionary-like object that contains information about the point that was clicked on
    :param figure: Figure object that is being updated
    :return: Updated figure object
    """
    if clickData:
        # Get the snp chr and pos and open the dbSNP page
        name = clickData['points'][0]['text']
        if '(' in name:
            snp = name.split('(')[1].rstrip(')')  # Extract the position information from the name string
        else:
            snp = name
        chromosome = snp.split(':')[0]
        position = snp.split(':')[1]

        # Open website
        #open_dbsnp_page(chromosome, position)
        open_ucsc_page(chromosome, position)

        point_index = clickData['points'][0]['pointIndex']
        colors = ['#ff0000' if i == point_index else 'blue' for i in range(len(df['MAF_to_plot']))]
        size = [10 if i == point_index else 15 for i in range(len(df['MAF_to_plot']))]
        figure['data'][0]['marker']['color'] = colors
        figure['data'][0]['marker']['size'] = size
    return figure


# Set up the layout
fig.update_layout(
    title='β vs. MAF Scatter Plot',
    xaxis=dict(title='MAF'),
    yaxis=dict(title='β')
)
app.layout = html.Div([
    dcc.Graph(id='scatter', figure=fig)
])

if __name__ == '__main__':
    app.run_server(debug=True)
