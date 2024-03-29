import os
import sys
import argparse
import pandas as pd
import dash
from dash import dcc, html
import plotly.graph_objs as go
import webbrowser

"""
This script creates an interactive plotly scatter plot from a BED-like tab-separated text file containing information 
about the MAFs, beta values and position of SNPs. It can also take the p-value and assigned gene for each variant.

The script visualizes the variants in a scatter plot, where the x-axis represents the MAFs and the y-axis represents the
beta values. The user can hover over any of the points to see more information about the variant. Clicking on a point
will open the USCS browser (https://genome-euro.ucsc.edu/cgi-bin/hgGateway) at the position of the variant. 
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
parser.add_argument('-pval', '--p_value_column', required=False,
                    help='The name of the column containing the p-values of the SNPs. Only used if specified.',
                    default=None)
parser.add_argument('-gene', '--gene_column', required=False,
                    help='The name of the column containing the gene names. Only used if specified.',
                    default=None)
parser.add_argument('-trait', '--trait_column', required=False,
                    help='The name of the trait. If none is specified, the default is "trait".',
                    default=None)
args = parser.parse_args()

# Check that the input file exists
if not os.path.isfile(args.input_file):
    raise ValueError('The input file does not exist')

# Read the input file
df = pd.read_csv(args.input_file, sep='\t')


def check_df_columns() -> None:
    """
    Check that the input file contains the needed columns, either specified by the user or the default column names.
    Note that 'name_column' is not required, as it can be created from the 'chromosome' and 'position' columns.

    :return:
    """
    for argument in [args.maf_column, args.beta_column, args.chromosome_column, args.position_column]:
        if argument not in df.columns:
            raise ValueError('The input file does not contain the column "{}"'.format(argument))
    if args.p_value_column and args.p_value_column not in df.columns:
        raise ValueError('The input file does not contain the column "{}"'.format(args.p_value_column))
    if args.gene_column and args.gene_column not in df.columns:
        raise ValueError('The input file does not contain the column "{}"'.format(args.gene_column))
    if args.name_column and args.name_column not in df.columns:
        raise ValueError('The input file does not contain the column "{}"'.format(args.name_column))


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
if args.trait_column:
    df['trait'] = df[args.trait_column].copy()
else:
    df['trait'] = None


def check_maf_percentage_or_fraction(in_df):
    """
    Crudely check if the MAF column is in percentages or fractions (e.g. 0.5 vs 50).
    If it is in percentages, convert it to fractions.
    This is a crude check that uses the ranges of the values to try to figure out what data we are dealing with.
        If the values range between 0 and 1, it might be fractions.
        If the values range between 0 and 100, it might be percentages.
    """
    # TODO


# For those entries where MAF > 0.5, flip it, and make change the sign of the beta
# TODO: be mindful that data can come both in percentages and fractions!!! (e.g. 0.5 vs 50). This should be harmonized!
mask = df[args.maf_column] > 50  # This only affects variants with MAF > 0.5
df.loc[mask, 'beta_to_plot'] = -df.loc[mask, 'beta_to_plot']
df.loc[mask, 'MAF_to_plot'] = 100 - df.loc[mask, 'MAF_to_plot']


# Keep only the variants with MAF < 2
# df = df[df['MAF_to_plot'] < 2]

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
    customdata=df[['p_value', 'gene', 'trait']],
)

# Create the figure object
fig = go.Figure(scatter)

# Define the layout hovertemplate based on the columns that were specified
# TODO: either do proper combinatorics or establish some sort of hierarchy for optional columns. In most cases all will be used anyway.
if args.p_value_column and args.gene_column and args.trait_column:
    hovertemplate_string = 'variant: %{text}<br>MAF: %{x}<br>β: %{y}<br>p-value: %{customdata[0]}<br>gene: %{customdata[1]}<br>trait: %{customdata[2]}'
elif args.p_value_column and args.gene_column:
    hovertemplate_string = 'variant: %{text}<br>MAF: %{x}<br>β: %{y}<br>p-value: %{customdata[0]}<br>gene: %{customdata[1]}'
elif args.p_value_column:
    hovertemplate_string = 'variant: %{text}<br>MAF: %{x}<br>β: %{y}<br>p-value: %{customdata[0]}'
elif args.gene_column:
    hovertemplate_string = 'variant: %{text}<br>MAF: %{x}<br>β: %{y}<br>gene: %{customdata[1]}'
else:
    hovertemplate_string = 'variant: %{text}<br>MAF: %{x}<br>β: %{y}'
fig.update_traces(hovertemplate=hovertemplate_string)


# Define callback to update scatter plot based on MAF threshold and handle click events
@app.callback(
    dash.dependencies.Output('scatter', 'figure'),
    [dash.dependencies.Input('maf-threshold-slider', 'value'),
     dash.dependencies.Input('scatter', 'clickData')],
    [dash.dependencies.State('scatter', 'figure')]
)
def update_scatterplot_and_point(maf_threshold, clickData, figure):
    """
    Update the scatter plot based on the MAF threshold and handle click events to open the desired url.
    """
    # Update the DataFrame based on the MAF threshold
    filtered_df = df.loc[df['MAF_to_plot'] < maf_threshold]
    # Recover the dimensions of the original DataFrame, just with None values where it was filtered out
    filtered_df = filtered_df.reindex(df.index)

    # Make a copy of the original figure to preserve the styling
    updated_scatter = figure.copy()

    # Update the scatter plot with the filtered DataFrame
    updated_scatter['data'][0]['x'] = filtered_df['MAF_to_plot']
    updated_scatter['data'][0]['y'] = filtered_df['beta_to_plot']
    updated_scatter['data'][0]['text'] = filtered_df['name']
    updated_scatter['data'][0]['customdata'][0] = filtered_df['p_value']
    updated_scatter['data'][0]['customdata'][1] = filtered_df['gene']
    updated_scatter['data'][0]['customdata'][2] = filtered_df['trait']

    # Handle click event
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
        # open_dbsnp_page(chromosome, position)
        open_ucsc_page(chromosome, position)

        point_index = clickData['points'][0]['pointIndex']
        colors = ['#ff0000' if i == point_index else 'blue' for i in range(len(df['MAF_to_plot']))]
        size = [10 if i == point_index else 15 for i in range(len(df['MAF_to_plot']))]
        updated_scatter['data'][0]['marker']['color'] = colors
        updated_scatter['data'][0]['marker']['size'] = size

    return updated_scatter


# Define callback to update click state
@app.callback(
    dash.dependencies.Output('scatter', 'clickData'),
    [dash.dependencies.Input('maf-threshold-slider', 'value')]
)
def reset_clickData(slider_value):
    """
    When the slider is moved, reset the clickData to None
    """
    return None


# Set up the layout
fig.update_layout(
    title='β vs. MAF Scatter Plot',
    xaxis=dict(title='MAF'),
    yaxis=dict(title='β')
)
app.layout = html.Div([
    dcc.Graph(id='scatter', figure=fig),
    dcc.Store(id='click-state', data={'clicked': False}),
    html.Div([
        dcc.Slider(
            id='maf-threshold-slider',
            min=0,
            max=50,
            step=1,
            marks={i: str(i) for i in range(50)},
            value=50,
            tooltip={'placement': 'bottom', 'always_visible': True},
        ),
    ], style={'margin': '300px'}),
])

if __name__ == '__main__':
    app.run(debug=True)

