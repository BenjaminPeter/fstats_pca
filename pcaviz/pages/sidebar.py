import dash
from dash import html
import dash_bootstrap_components as dbc


def nav():
    nav = dbc.Nav(
    [
        dbc.NavLink([html.Div(page['name'], className='ms-2'),],
                    href=page['path'],
                    active='exact',)
        for page in dash.page_registry.values()
    ],
        pills=True,
        className = 'bg-light'
    )
    return nav


def sidebar():
    sb =  html.Div(
        dbc.Nav(
            [
                dbc.NavLink(
                    [
                        html.Div(page["name"], className="ms-2"),
                    ],
                    href=page["path"],
                    active="exact",
                )
                for page in dash.page_registry.values()
            ],
            vertical=True,
            pills=True,
            className="bg-light",
        )
    )

    return sb

