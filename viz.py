from dash import Dash, html, dcc, Input, Output, ctx
from dash.exceptions import PreventUpdate
import pandas as pd
import plotly.express as px

external_stylesheets = []

app = Dash(__name__, external_stylesheets=external_stylesheets)

df2 = pd.read_csv('P0.csv')
pcs = [c for c in df2.columns if c.startswith("PC")]
meta_cols = [c for c in df2.columns if not c.startswith("PC")]
pops = pd.unique(df2['pop'])
df2.head()


app.layout = html.Div([
    html.Div([

        html.Div([
            dcc.Dropdown(
                pcs,
                pcs[0],
                id='crossfilter-xaxis-column',
            )
        ],
        style={'width': '49%', 'display': 'inline-block'}),

        html.Div([
            dcc.Dropdown(
                pcs,
                pcs[1],
                id='crossfilter-yaxis-column'
            ),
        ], style={'width': '49%', 'float': 'right', 'display': 'inline-block'}),
        html.Div([
            dcc.Dropdown(
                pops,
                pops[0],
                id='selector-f2-pop1'
            ),
        ], style={'width': '49%', 'display': 'inline-block'}),
        html.Div([
            dcc.Dropdown(
                pops,
                pops[1],
                id='selector-f2-pop2'
            ),
        ], style={'width': '49%', 'float': 'right', 'display': 'inline-block'})
    ]),

    html.Div([
        dcc.Graph(
            id='pcplot',
            hoverData={'points': [{'customdata': 'PC3'}]}
        )
    ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),
    html.Div([
        #dcc.Graph(id='x-time-series'),
        #dcc.Graph(id='y-time-series'),
    ], style={'display': 'inline-block', 'width': '49%'}),

    #html.Div(dcc.Slider(
    #    df2.columns[0],
    #    df2.columns[10],
    #    step=None,
    #    id='crossfilter-year--slider',
    #    value=df2.columns[0],
    #    marks=df2.columns[:10]
    #), style={'width': '49%', 'padding': '0px 20px 20px 20px'})
])


@app.callback(
    Output('pcplot', 'figure'),
    Input('crossfilter-xaxis-column', 'value'),
    Input('crossfilter-yaxis-column', 'value'))
#    Input('crossfilter-year--slider', 'value'))
def update_graph(xaxis_column_name, yaxis_column_name):

    fig = px.scatter(df2,
        x=xaxis_column_name,
        y=yaxis_column_name,
        color = 'pop',
        hover_name='sample',
        hover_data = meta_cols + [xaxis_column_name, yaxis_column_name]
            )

    #fig.update_traces(customdata=dff[dff['Indicator Name'] == yaxis_column_name]['Country Name'])

    fig.update_xaxes(title=xaxis_column_name, type='linear')
    fig.update_yaxes(title=yaxis_column_name, type='linear')

    fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')

    return fig


def create_time_series(dff, title):
    return None

    fig = px.scatter(dff, x='Year', y='Value')

    fig.update_traces(mode='lines+markers')

    fig.update_xaxes(showgrid=False)

    fig.update_yaxes(type='linear')

    fig.add_annotation(x=0, y=0.85, xanchor='left', yanchor='bottom',
                       xref='paper', yref='paper', showarrow=False, align='left',
                       text=title)

    fig.update_layout(height=225, margin={'l': 20, 'b': 30, 'r': 10, 't': 10})
    return None

    return fig

@app.callback(
    Output('selector-f2-pop1', 'value'),
    Input('pcplot', 'clickData')
)
def update_f2_dropdown1(click_data):
    if not click_data:
        raise PreventUpdate
    breakpoint()
    p = click_data['points']
    clicked_sample = p[0]['hovertext']
    clicked_sample_data = df2[df2['sample'] == clicked_sample]
    clicked_pop = clicked_sample_data['pop']
    return clicked_pop.iloc[0]


if False:
    @app.callback(
        Output('x-time-series', 'figure'),
        Input('', 'hoverData'),
        Input('crossfilter-xaxis-column', 'value'),
    )
    def update_y_timeseries(hoverData, xaxis_column_name):
        country_name = hoverData['points'][0]['customdata']
        dff = df[df['Country Name'] == country_name]
        dff = dff[dff['Indicator Name'] == xaxis_column_name]
        title = '<b>{}</b><br>{}'.format(country_name, xaxis_column_name)
        return create_time_series(dff, title)


    @app.callback(
        Output('y-time-series', 'figure'),
        Input('', 'hoverData'),
        Input('crossfilter-yaxis-column', 'value'),
    )
    def update_x_timeseries(hoverData, yaxis_column_name, axis_type):
        dff = df[df['Country Name'] == hoverData['points'][0]['customdata']]
        dff = dff[dff['Indicator Name'] == yaxis_column_name]
        return create_time_series(dff, yaxis_column_name)

if __name__ == '__main__':
    app.run_server(debug=True)

