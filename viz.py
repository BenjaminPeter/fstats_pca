from dash import Dash, html, dcc, Input, Output, ctx
from dash.exceptions import PreventUpdate
import pandas as pd
import plotly.express as px
import numpy as np

external_stylesheets = []

app = Dash(__name__, external_stylesheets=external_stylesheets)

df2 = pd.read_csv("P0.csv")
pcs = [c for c in df2.columns if c.startswith("PC")]
df = pd.melt(id_vars=["sample", "pop"], value_vars=pcs, frame=df2, var_name="pc")
meta_cols = [c for c in df2.columns if not c.startswith("PC")]
pops = pd.unique(df2["pop"])


app.layout = html.Div(
    [
        html.Div( #top bar with PCs and population selectors
            [
                html.Div(
                    [
                        dcc.Dropdown(pops, pops[0], id="selector-f2-pop1"),
                    ],
                    style={
                        "width": "24%",
                        "float": "left",
                        "display": "inline-block",
                    },
                ),
                html.Div(
                    [
                        dcc.Dropdown(
                            pops, [pops[1], pops[2]], multi=True, id="selector-f2-pop2"
                        ),
                    ],
                    style={
                        "width": "24%",
                        "float": "left",
                        "display": "inline-block",
                    },
                ),
                html.Div(
                    [
                        dcc.Dropdown(
                            pops, [pops[3], pops[4]], multi=True, id="selector-f2-pop3"
                        ),
                    ],
                    style={
                        "width": "24%",
                        "float": "left",
                        "display": "inline-block",
                    },
                ),
                html.Div(
                    [
                        dcc.Dropdown(
                            pops, [pops[-1]], multi=True, id="selector-f2-pop4"
                        ),
                    ],
                    style={
                        "width": "24%",
                        "float": "left",
                        "display": "inline-block",
                    },
                ),
                html.Div(
                    dcc.Slider(
                        1,
                        len(pcs),
                        1,
                        id="slider-n-pcs",
                        value=10,
                        marks=dict(tpl for tpl in enumerate(pcs) if tpl[0] % 5 == 4),
                    ),
                    style={
                        "width": "49%",
                        "float": "right",
                        "width": "99%",
                        "padding": "40px 20px 20px 20px",
                    },
                ),
            ],
            style={
                "width": "100%",
                "float": "right",
                "border": "1px red solid",
            },
        ),
        html.Div( #main plot div
            [
                dcc.Graph(id="pcplot", hoverData={"points": [{"customdata": "PC3"}]}),
                html.Div(
                    [
                        "X-Axis",
                        dcc.Dropdown(
                            pcs,
                            pcs[0],
                            clearable=False,
                            id="xaxis-pca",
                        ),
                    ],
                    style={"width": "19%", "float": "left", "display": "inline-block"},
                ),
                html.Div(
                    [
                        "Y-Axis",
                        dcc.Dropdown(pcs, pcs[1], clearable=False, id="yaxis-pca"),
                    ],
                    style={"width": "19%", "float": "left", "display": "inline-block"},
                ),
                html.Div(
                    [
                        "Clicking will affect...",
                        dcc.RadioItems(
                            ["Pos1", "Pos2", "Pos3", "Pos4"],
                            "Pos1",
                            inline=True,
                            id="radio-click",
                        ),
                    ],
                    style={
                        "width": "59%",
                        "float": "bottom left",
                        "display": "inline-block",
                    },
                ),
            ],
            style={
                "width": "59%",
                "border": "1px black solid",
                "display": "inline-block",
                "padding": "0 20",
            },
        ),
        html.Div( #left plots div
            [
                dcc.Graph(id="f2-plot1"),
                dcc.Graph(id="f2-plot2"),
            ],
            style={
                "display": "inline-block",
                "width": "39%",
                "border": "1px black solid",
            },
        ),
    ]
)


@app.callback(
    Output("pcplot", "figure"), Input("xaxis-pca", "value"), Input("yaxis-pca", "value")
)
#    Input('crossfilter-year--slider', 'value'))
def update_pca_biplot(xaxis_column_name, yaxis_column_name):
    fig = px.scatter(
        df2,
        x=xaxis_column_name,
        y=yaxis_column_name,
        color="pop",
        hover_name="sample",
        hover_data=meta_cols + [xaxis_column_name, yaxis_column_name],
    )

    fig.update_xaxes(title=xaxis_column_name, type="linear")
    fig.update_yaxes(
        title=yaxis_column_name, type="linear", scaleanchor="x", scaleratio=1
    )

    fig.update_layout(margin={"l": 40, "b": 40, "t": 10, "r": 0}, hovermode="closest")

    return fig


def create_f2_plot_by_pcs(f2_p, n_pcs):
    """returns a plot that has PCs on X and f2-values on Y"""

    fig = px.bar(f2_p, x="value", color="pops", facet_col="pops", error_x="error")

    # fig.update_traces(mode='lines+markers')

    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(categoryorder="array", categoryarray=pcs[:n_pcs])

    fig.update_layout(height=225, margin={"l": 20, "b": 30, "r": 10, "t": 10})

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

    fig.update_layout(height=225, margin={"l": 20, "b": 30, "r": 10, "t": 10})

    return fig


@app.callback(Output("selector-f2-pop1", "value"), Input("pcplot", "clickData"))
def update_f2_dropdown1(click_data):
    """basic update to have clicked plot be first ind"""
    if not click_data:
        raise PreventUpdate
    p = click_data["points"]
    clicked_sample = p[0]["hovertext"]
    clicked_sample_data = df2[df2["sample"] == clicked_sample]
    clicked_pop = clicked_sample_data["pop"]
    return clicked_pop.iloc[0]


def make_fstat_table(pop1, pop2, n_pcs):
    if type(pop1) is str:
        pop1 = [pop1]
    if type(pop2) is str:
        pop2 = [pop2]

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


@app.callback(
    Output("f2-plot1", "figure"),
    Input("selector-f2-pop1", "value"),
    Input("selector-f2-pop2", "value"),
    Input("slider-n-pcs", "value"),
)
def update_f2_plot_by_pc(pop1, pop2, n_pcs):
    dfp = make_fstat_table(pop1, pop2, n_pcs)
    return create_f2_plot_by_pcs(dfp, n_pcs)


@app.callback(
    Output("f2-plot2", "figure"),
    Input("selector-f2-pop1", "value"),
    Input("selector-f2-pop2", "value"),
    Input("slider-n-pcs", "value"),
)
def update_f2_plot_series(pop1, pop2, n_pcs):
    dfp = make_fstat_table(pop1, pop2, n_pcs)
    return create_f2_plot_series(dfp, n_pcs)


if __name__ == "__main__":
    app.run_server(debug=True)
