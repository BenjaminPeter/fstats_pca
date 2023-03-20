import dash
from dash import Dash, html, dcc, Input, Output, State, callback
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import logging





logger = logging.getLogger("")
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler("spam.log")
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)

def norm2(v):
    return np.sum(v ** 2)

def proj(A, v):
    Av = np.dot(A, v)
    norm = norm2(v)
    proj = np.outer(Av, v) / norm
    return proj



def load_data():
    df2 = pd.read_csv("P0.csv")
    df_pop = df2.groupby("pop").mean().reset_index()
    pcs = [c for c in df2.columns if c.startswith("PC")]
    df = pd.melt(id_vars=["sample", "pop"], value_vars=pcs, frame=df2, var_name="pc")
    meta_cols = [c for c in df2.columns if not c.startswith("PC")]
    pops = pd.unique(df2["pop"])
    inds = pd.unique(df2["sample"])

    ev = df2[pcs].apply(lambda x: sum(x**2))
    pcte = ev / sum(ev)
    return df, df2, pcs, pops, inds, meta_cols, pcte, df_pop


df, df2, pcs, pops, inds, meta_cols, pcte, df_pop = load_data()
dfi = df_pop.set_index('pop')
dfi = dfi.iloc[:10, :2]
dfv = dfi.values
proj_axis = np.diff(dfi.loc[['Pop 0', 'Pop 1']],axis=0).T
proj0 = dfi @ proj_axis / norm2(proj_axis)
r0 = dfi - (proj0 @ proj_axis.T).values



cols = px.colors.qualitative.Alphabet
col_dict = dict( (pop, cols[i % len(cols)]) for i, pop in enumerate(pops))


dash.register_page(__name__, name="Pca Viz1", top_nav=True)

#breakpoint()
def proj(dfi, proj_pops=['Pop 0', 'Pop 1']):
    proj_axis = np.diff(dfi.loc[proj_pops],axis=0)
    pass



def layout():
    layout = html.Div(
        [
            html.Div(  # top bar with PCs and population selectors
                layout_selector(),
                style={
                    "width": "100%",
                    "float": "right",
                    "border": "1px red solid",
                },
            ),
            html.Div(
                layout_pca_biplot(),
                style={
                    "width": "59%",
                    "border": "1px black solid",
                    "display": "inline-block",
                    "padding": "0 20",
                    "verticalAlign": "top",
                },
            ),
            html.Div(  # left plots div
                layout_right_plots(),
                style={
                    "display": "inline-block",
                    "width": "39%",
                    "border": "1px black solid",
                },
            ),
        ]
    )
    return layout

def layout_selector():
    return [
        html.Div(  # pop1
            [
                dcc.Dropdown(pops, [pops[0]], multi=True, id="selector-f2-pop1"),
                dcc.Dropdown(inds, [inds[0]], multi=True, id="selector-f2-ind1"),
            ],
            style={
                "width": "19%",
                "float": "left",
                "display": "inline-block",
            },
        ),
        html.Div(
            [
                dcc.Dropdown(
                    pops, [pops[1], pops[2]], multi=True, id="selector-f2-pop2"
                ),
                dcc.Dropdown(inds, [inds[0]], multi=True, id="selector-f2-ind2"),
            ],
            style={
                "width": "19%",
                "float": "left",
                "display": "inline-block",
            },
        ),
        html.Div(
            [
                dcc.Dropdown(
                    pops, [pops[3], pops[4]], multi=True, id="selector-f2-pop3"
                ),
                dcc.Dropdown(inds, [inds[0]], multi=True, id="selector-f2-ind3"),
            ],
            style={
                "width": "19%",
                "float": "left",
                "display": "inline-block",
            },
        ),
        html.Div(
            [
                dcc.Dropdown(pops, [pops[-1]], multi=True, id="selector-f2-pop4"),
                dcc.Dropdown(inds, [inds[0]], multi=True, id="selector-f2-ind4"),
            ],
            style={
                "width": "19%",
                "float": "left",
                "display": "inline-block",
            },
        ),
        html.Div([
            dbc.Button("Set #PCS", id='button-n-slider', color='primary',
                       n_clicks=0)
        ],
            style={
                "width": "19%",
                "float": "left",
                "display": "inline-block",
            },
         ),
        html.Div(
            dbc.Collapse([
            dcc.Slider(
                0,
                len(pcs),
                1,
                id="slider-n-pcs",
                value=9,
                marks=dict(tpl for tpl in enumerate(pcs) if tpl[0] % 5 == 4),
            ),
            ]
                         ,id = 'collapse-n-pcs',
                         is_open=False,
                         ),
            style={
                "width": "49%",
                "float": "right",
                "width": "99%",
                "padding": "40px 20px 20px 20px",
            },
        ),
    ]


def layout_pca_biplot():
    return [
        dcc.Graph(
            id="pcplot",
            responsive=True,
            style={"height": "70vh"},
            hoverData={"points": [{"customdata": "PC3"}]},
        ),
        html.Div(
            [
                dbc.Label("X-Axis"),
                dcc.Dropdown(
                    pcs,
                    pcs[0],
                    clearable=False,
                    id="xaxis-pca",
                ),
                dbc.Label("Y-Axis"),
                dcc.Dropdown(pcs, pcs[1], clearable=False, id="yaxis-pca"),
            ],
            style={"width": "19%", "float": "left", "display": "inline-block"},
        ),
        html.Div(
            [
                dbc.Label("Clicking will affect..."),
                dbc.RadioItems(
                    ["Pos1", "Pos2", "Pos3", "Pos4"],
                    "Pos1",
                    inline=False,
                    id="radio-click-target",
                ),
            ],
            style={
                "width": "19%",
                "float": "bottom left",
                "display": "inline-block",
                "padding": "0px 20px 20px 20px",
            },
        ),
        html.Div(
            [
                dbc.Label("Visuals to display..."),
                dbc.Checklist(
                    options=[
                        {"label": "Show Individual-F2", "value": 'f2ind',
                         'disabled': True},
                        {"label": "Show Pop-F2", "value": 'f2pop', "disabled":
                         False},
                        {"label": "Show Ind PC", "value": 'pcind', "disabled":
                         False},
                        {"label": "Show Pop PC", "value": 'pcpop', "disabled": False},
                        {"label": "Show Admixture F3 circle", "value":
                         'af3circle', "disabled": False},
                    ],
                    inline = True,
                    value=['pcpop', 'f2pop', 'af3circle'],
                    id="checklist-input",
                )
            ],
            style={
                "width": "29%",
                "float": "bottom left",
                "display": "inline-block",
                "padding": "0px 20px 20px 20px",
            },
        ),
    ]

def layout_f2_tab():
    return [
        dcc.Graph(id="f2-plot1"),
        dcc.Graph(id="f2-plot2"),
    ]

def layout_af3_tab():
    return [
        dcc.Graph(id="f3-plot1"),
        dcc.Graph(id="f3-plot2"),
    ]

def layout_right_plots():
    tabs = html.Div([
        dbc.Tabs([
        dbc.Tab(label="F2", tab_id='tab-f2'),
        dbc.Tab(label="Admixture-F3", tab_id='tab-af3'),

    ], active_tab='tab-f2', id='tabs-right'),
        html.Div(id='right-content')
    ])

    return tabs

@callback(Output("right-content", "children"), [Input("tabs-right", "active_tab")])
def switch_tab(at):
    if at == "tab-f2":
        return layout_f2_tab()
    elif at == "tab-af3":
        return layout_af3_tab()
    return html.P("This shouldn't ever be displayed...")

@callback(
    Output("collapse-n-pcs", "is_open"),
    [Input("button-n-slider", "n_clicks")],
    [State("collapse-n-pcs", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@callback(
    Output("pcplot", "figure"),
    Input("xaxis-pca", "value"),
    Input("yaxis-pca", "value"),
    Input("selector-f2-pop1", "value"),
    Input("selector-f2-pop2", "value"),
    Input("selector-f2-pop3", "value"),
    Input("selector-f2-pop4", "value"),
    Input("checklist-input", "value"),
)
#    Input('crossfilter-year--slider', 'value'))
def update_pca_biplot(pcx, pcy, pop1, pop2, pop3,
                      pop4, opt):
    fig = go.Figure()
    if 'pcind' in opt:
        s1 = px.scatter(
            df2,
            x=pcx,
            y=pcy,
            color="pop",
            color_discrete_map = col_dict,
            hover_name="sample",
            hover_data=meta_cols + [pcx, pcy],
        )
        s1.update_traces(marker=dict(symbol='circle-open'))
        for pts in s1.data:
            fig.add_trace(pts)
    if 'pcpop' in opt:
        s2 = px.scatter(
            df_pop,
            x=pcx,
            y=pcy,
            color="pop",
            color_discrete_map = col_dict,
            hover_name="pop",
            hover_data= [pcx, pcy],
        )
        s2.update_traces(marker=dict(symbol='square', size=10))
        
        for pts in s2.data:
            pts.showlegend=False
            fig.add_trace(pts)
    if 'f2pop' in opt:
        for p1 in pop1:
            for p2 in pop2:
                xpts, ypts = [], []
                xpts.append(df_pop[df_pop['pop'] ==
                                   p1][pcx].iloc[0])
                xpts.append(df_pop[df_pop['pop'] ==
                                   p2][pcx].iloc[0])
                ypts.append(df_pop[df_pop['pop'] ==
                                   p1][pcy].iloc[0])
                ypts.append(df_pop[df_pop['pop'] ==
                                   p2][pcy].iloc[0])
                s3 = go.Scatter(x=xpts, y=ypts, mode = 'lines', 
                                line_dash='dash',
                                line_color='black',
                                line_width=1,
                                showlegend=False,
                                opacity=.5)
                fig.add_trace(s3)
    if 'af3circle' in opt:
        dfp = df_pop.set_index('pop')[[pcx, pcy]]
        for p2 in pop1:
            for p3 in pop2:
                center = dfp.loc[[p2,p3]].mean()
                k = df_pop.set_index('pop').loc[[p2,p3]]
                r1 = np.sqrt(np.sum((k.iloc[0] - k.iloc[1]) ** 2)) /2
                x0, x1  = center.iloc[0] - r1, center.iloc[0] + r1
                y0, y1  = center.iloc[1] - r1, center.iloc[1] + r1
                fig.add_shape(type="circle",
                    xref="x", yref="y",
                              x0=x0, y0=y0, x1=x1, y1=y1,
    line_color="LightSeaGreen",
)








    xlab = f"{pcx} ({100*pcte.loc[pcx]:2.2f}%)"
    ylab = f"{pcy} ({100*pcte.loc[pcy]:2.2f}%)"
    fig.update_xaxes(title=xlab)
    fig.update_yaxes(title=ylab, scaleanchor="x", scaleratio=1)

    f2_viz = px.line(x=[0, 100, 200], y=[0, -100, 200])
    fig.update_layout(
        margin={"l": 40, "b": 40, "t": 10, "r": 10}, autosize=True, hovermode="closest"
    )

    return fig


def create_f2_plot_by_pcs(f2_p, n_pcs):
    """returns a plot that has PCs on X and f2-values on Y"""

    fig = px.bar(f2_p, x="value", color="pops", facet_col="pop_x", 
                 facet_row='pop_y',
                 error_x="error")

    # fig.update_traces(mode='lines+markers')

    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(categoryorder="array", categoryarray=pcs[:n_pcs])

    fig.update_layout(height=300, 
                      showlegend=False,
                      margin={"l": 20, "b": 30, "r": 10, "t": 30})

    return fig


def create_f2_plot_series(f2_p, n_pcs):
    """returns a plot that has f2 on X and pops on y"""

    df = f2_p.groupby("pops").sum()

    fig = px.bar(df, x="value", error_x="error")

    # fig.update_traces(mode='lines+markers')

    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(categoryorder="array", categoryarray=pcs[:n_pcs])

    # fig.add_annotation(x=0, y=0.85, xanchor='left', yanchor='bottom',
    #                   xref='paper', yref='paper', showarrow=False, align='left',
    #                   text=title)

    fig.update_layout(height=300, margin={"l": 20, "b": 30, "r": 10, "t": 10})

    return fig


@callback(
    inputs=[
        Input("pcplot", "clickData"),
        State("radio-click-target", "value"),
        State("selector-f2-pop1", "value"),
        State("selector-f2-pop2", "value"),
        State("selector-f2-pop3", "value"),
        State("selector-f2-pop4", "value"),
    ],
    output=[
        Output("selector-f2-pop1", "value"),
        Output("selector-f2-pop2", "value"),
        Output("selector-f2-pop3", "value"),
        Output("selector-f2-pop4", "value"),
        Output("pcplot", "clickData"),  # required for repeated clicks
    ],
    prevent_initial_call=True,
)
def update_f2_dropdown1(click_data, click_target, *selectors):
    """basic update to have clicked plot be first ind"""
    logger.debug(click_data)
    logger.debug(selectors)
    if not click_data:
        raise PreventUpdate
    p = click_data["points"]

    clicked_name = p[0]["hovertext"]
    if clicked_name in pops:
        clicked_pop = clicked_name
    else:
        clicked_sample_data = df2[df2["sample"] == clicked_name]
        clicked_pop = clicked_sample_data["pop"].iloc[0]
    target = int(click_target[3]) - 1
    if clicked_pop in selectors[target]:
        selectors[target].remove(clicked_pop)
    else:
        selectors[target].append(clicked_pop)
    logger.debug(selectors)
    return selectors + (None,)


def make_f2_table(pop1, pop2, n_pcs):
    df1p = df[(df["pop"].isin(pop1)) & (df["pc"].isin(pcs[:n_pcs]))]
    x1 = df1p.groupby(["pop", "pc"]).mean().reset_index("pop")

    df2p = df[(df["pop"].isin(pop2)) & (df["pc"].isin(pcs[:n_pcs]))]
    x2 = df2p.groupby(["pop", "pc"]).mean().reset_index("pop")

    dff = x1.merge(x2, on="pc")
    dfp = pd.DataFrame()
    dfp["value"] = (dff.value_x - dff.value_y) ** 2
    dfp["pop_x"] = dff.pop_x
    dfp["pop_y"] = dff.pop_y
    dfp["pops"] = "f2(" + dfp.pop_x + ", " + dfp.pop_y + ")"
    dfp["error"] = dfp["value"] / 10  # dummy error

    return dfp


@callback(
    Output("f2-plot1", "figure"),
    Input("selector-f2-pop1", "value"),
    Input("selector-f2-pop2", "value"),
    Input("slider-n-pcs", "value"),
)
def update_f2_plot_by_pc(pop1, pop2, n_pcs):
    dfp = make_f2_table(pop1, pop2, n_pcs+1)
    return create_f2_plot_by_pcs(dfp, n_pcs+1)


@callback(
    Output("f2-plot2", "figure"),
    Input("selector-f2-pop1", "value"),
    Input("selector-f2-pop2", "value"),
    Input("slider-n-pcs", "value"),
)
def update_f2_plot_series(pop1, pop2, n_pcs):
    dfp = make_f2_table(pop1, pop2, n_pcs+1)
    return create_f2_plot_series(dfp, n_pcs+1)


if __name__ == "__main__":
    app.run_server(debug=True)
