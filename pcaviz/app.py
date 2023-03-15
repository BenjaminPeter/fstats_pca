from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
import dash
from pages.sidebar import sidebar, nav

app = Dash(__name__, use_pages=True, 
           external_stylesheets=[dbc.themes.BOOTSTRAP],
           suppress_callback_exceptions=True)

app.layout = html.Div([
    nav(),
	dash.page_container
])

if __name__ == '__main__':
	app.run_server(debug=True)
