import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
from .sidebar import sidebar, nav

dash.register_page(__name__, path='/', top_nav=True)





xxlayout0 = html.Div(children=[
    html.H1(children='This is our Home page'),

    html.Div(children='''
        This is our Home page content.
    '''),

])

def layout():
    return 'HOME' #[nav()]#, dbc.Col(sidebar(), width=2), dbc.Col(layout0, width=10)]
