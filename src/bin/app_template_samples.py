import dash
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests as req

from dash import dcc, html
from dash.dependencies import Input, Output
from dotenv import load_dotenv

app = dash.Dash(__name__)
ref = "{replace_ref}"
sitename = "{replace_sitename}"
cov_df = pd.read_csv(f"{sitename}.{ref}__all.finalCov.tsv", sep="\t")
region_df = pd.read_csv("{regionsStandardizedTsv}", sep="\t")
region_df = (
    region_df[region_df.chrom == ref]
    .sort_values(by="st")[["chrom", "st", "end", "region"]]
    .drop_duplicates(subset=["chrom", "st", "end"])
)
region_ls = ["All"] + list(set(region_df["region"].values))
region_y_last_end = -1
region_y_last = 1
region_y_base = 1
display_y = []
tmp_end_ls = [region_y_last_end]


def pass_ends(st, end_ls):
    for end in end_ls:
        if st - 10 < end:
            return False
    return True


for _, row in region_df.iterrows():
    end = row.end
    if pass_ends(row.st, tmp_end_ls):
        # reset y
        use_y = region_y_base
    else:
        use_y = region_y_last + 1
    region_y_last = use_y
    region_y_last_end = end
    tmp_end_ls.append(region_y_last_end)
    display_y.append(use_y)

region_df["region_y"] = display_y
max_region_y = max(region_df["region_y"])

x_st = min(cov_df["pos"])
x_end = max(cov_df["pos"])


location_options = [{"label": r, "value": r} for r in region_ls]

app.layout = html.Div(
    [
        dcc.Store(id="zoom-st"),
        dcc.Store(id="zoom-end"),
        html.Div(
            [
                html.Div(
                    [
                        html.Img(
                            src="https://d33wubrfki0l68.cloudfront.net/1ac3f0e3753f18c7e2a8893957d1841fba1e3d08/48a33/wp-content/uploads/2018/10/rstudio-logo-flat.png",
                            style={"height": "60px", "width": "auto"},
                        )
                    ],
                    className="one-third column",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H3(
                                    "All samples",
                                    style={"margin-bottom": "0px"},
                                ),
                                html.H5(ref, style={"margin-top": "0px"}),
                            ]
                        )
                    ],
                    className="one-half column",
                    id="title",
                ),
            ],
            id="header",
            className="row flex-display",
            style={"margin-bottom": "25px"},
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.P("Select Region:"),
                                dcc.Dropdown(
                                    id="location", options=location_options, value="All"
                                ),
                            ],
                            style={"margin-top": "10"},
                        ),
                    ],
                    className="row",
                ),
                html.Div(
                    [
                        html.Div(
                            [dcc.Graph(id="coverage-plot")],
                            style={"margin-bottom": "1"},
                        ),
                        html.Div(
                            [dcc.Graph(id="region-map")], style={"margin-top": "1"}
                        ),
                    ]
                ),
            ],
            className="row",
        ),
    ]
)


@app.callback(Output("coverage-plot", "figure"), [Input("location", "value")])
def update_forecast_graph(value):
    if value == "All" or not value:
        subset_df = cov_df
    else:
        # allow multiple rows
        st = list(region_df[region_df["region"] == value]["st"])[0]
        end = list(region_df[region_df["region"] == value]["end"])[0]
        subset_df = cov_df[(cov_df.pos >= st) & (cov_df.pos <= end)]

    fig = px.line(subset_df, x="pos", y="depth", color="sample_name")
    fig.update_layout(showlegend=False)
    fig.update_layout(
        xaxis_title="Position",
        yaxis_title="Read count",
    )
    return fig


@app.callback(
    Output("region-map", "figure"),
    [Input("location", "value"), Input("coverage-plot", "relayoutData")],
)
def update_region_graph(value, relayoutData):
    init_region = False
    if dash.callback_context.triggered:
        for cb in dash.callback_context.triggered:
            if cb["prop_id"] == "location.value" and not value:
                init_region = True

    fig = go.Figure()

    x_trace, y_trace, label_trace = [], [], []
    for _, row in region_df.iterrows():
        x0 = row.st
        x1 = row.end
        y0 = row.region_y + 0.25
        y1 = row.region_y - 0.25
        x_trace.append(x0)
        y_trace.append(y0)
        label_trace.append(row.region)
        fig.add_trace(
            go.Scatter(
                x=[x0, x0, x1, x1],
                y=[y0, y1, y1, y0],
                fill="toself",
                hoveron="points+fills",
                name="",
                text=row.region,
                opacity=0,
                hoverinfo="text",
            )
        )

    fig.add_trace(
        go.Scatter(
            x=x_trace,
            y=y_trace,
            text=label_trace,
            mode="text",
        )
    )

    # Set axes properties
    has_coords = False
    if not init_region:
        if relayoutData:
            if "xaxis.range[0]" in relayoutData:
                st = relayoutData["xaxis.range[0]"]
                end = relayoutData["xaxis.range[1]"]
                has_coords = True
            elif "autosize" in relayoutData and not value:
                st = x_st
                end = s_end
                has_coords = True

    if not has_coords:
        # assign from zoom region
        if value == "All" or not value:
            st = x_st
            end = x_end
        else:
            st = list(region_df[region_df["region"] == value]["st"])[0]
            end = list(region_df[region_df["region"] == value]["end"])[0]

    fig.update_xaxes(range=[st, end], showgrid=False)
    fig.update_yaxes(range=[0, max_region_y + 0.5], visible=False, showgrid=False)

    # Add shapes
    for _, row in region_df.iterrows():
        x0 = row.st
        x1 = row.end
        y0 = row.region_y
        y1 = row.region_y
        fig.add_shape(
            type="rect",
            x0=x0,
            y0=y0,
            x1=x1,
            y1=y1,
            line=dict(color="RoyalBlue"),
        )

    fig.update_shapes(dict(xref="x", yref="y"))
    fig.update_layout(
        {
            "plot_bgcolor": "rgba(0, 0, 0, 0)",
            "paper_bgcolor": "rgba(0, 0, 0, 0)",
            "showlegend": False,
        }
    )
    return fig


if __name__ == "__main__":
    app.run_server(debug=True)
