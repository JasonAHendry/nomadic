import os
import logging
import datetime
import pandas as pd
import seaborn as sns

import plotly.graph_objects as go

from matplotlib.colors import rgb2hex

from dash import Dash, html, dcc
from dash.dependencies import Input, Output



def create_dashboard(
        expt_name: str,
        flagstats_csv: str,
        bedcov_csv: str
    ):
    """
    Create a dashboard monitoring the indicated files
    
    """

    # PARAMETEERS
    CATEGORIES = ["n_uniq_mapped", "n_chim_mapped", "n_mult_mapped", "n_unmapped"]
    COLS = [rgb2hex(c) for c in sns.color_palette("Blues_r") + [(0.75, 0.75, 0.75)]]
    COL_MAP = dict(zip(CATEGORIES, COLS))

    # Launch time
    t0 = datetime.datetime.now().replace(microsecond=0)

    # Create
    app = Dash(__name__, external_stylesheets=["assets/dashboard-style.css"])
    app_log = logging.getLogger("werkzeug")
    app_log.setLevel(logging.ERROR)

    top_row = html.Div(
        className='top-row',
        children=[
            html.Img(id="logo", src="assets/nomadic_logo-01.png"),
            html.Div(id="expt-summary")
        ])
    middle_row = html.Div(
        className="middle-row",
        children=[
            dcc.Graph(id='mapping-pie'),
            dcc.Graph(id='mapping-barplot')
        ]
    )
    bottom_row = html.Div(
        className="bottom-row",
        children=[
            dcc.Graph(id="bedcov-pie"),
            dcc.Graph(id="bedcov-stripplot")
        ]
    )

    layout = [top_row,
              html.Hr(),
              html.H2("Read Mapping Statistics"),
              middle_row,
              dcc.Checklist(id='mapping-checklist',
                            options=CATEGORIES,
                            value=CATEGORIES,
                            inline=True),
            html.Hr(),
            html.H2("Region Coverage Statistics"),
            dcc.Dropdown(id='bedcov-dropdown',
            options=["mean_cov", "n_reads", "cov_gr100_per"],
            value="n_reads"),
            bottom_row,
            html.Hr(),
            html.Br(),
            dcc.Interval(
                id='interval',
                interval=1_000,
                n_intervals=0
            )
    ]

    @app.callback(
        Output('bedcov-pie', 'figure'),
        Input('interval', 'n_intervals'),
        Input('bedcov-dropdown', 'value')
    )
    def update_bedcov_pie(_, dropdown_stat):
        """
        Update the pie chart

        """
        if not os.path.exists(bedcov_csv):
            return go.Figure()  # returing blank figure

        # Load data
        df = pd.read_csv(bedcov_csv)

        # Pie values
        pie_data = df.groupby("name", sort=True)[dropdown_stat].sum()
       # pie_values = df[dropdown_stat].sum().tolist()

        # Plot data
        plot_data = [
            go.Pie(labels=pie_data.index,
                   values=pie_data.values)
        ]
        fig = go.Figure(data=plot_data)

        return fig

    @app.callback(
        Output('bedcov-stripplot', 'figure'),
        Input('interval', 'n_intervals'),
        Input('bedcov-dropdown', 'value')
    )
    def update_stripplot(_, dropdown_stat):
        """
        Update the BED coverage stripplot
        
        """
        
        if not os.path.exists(bedcov_csv):
            return go.Figure()  # returing blank figure

        # Load data
        df = pd.read_csv(bedcov_csv)

        plot_data = [
            go.Scatter(
                x=tdf["barcode"],
                y=tdf[dropdown_stat],
                mode="markers",
                name=target
            )
            for target, tdf in df.groupby("name", sort=True)
        ]

        fig = go.Figure()
        for plot_trace in plot_data:
            fig.add_trace(plot_trace)

        fig.update_layout(
            yaxis_title=dropdown_stat,
            xaxis=dict(showline=True, linewidth=1, linecolor='black', mirror=True),
            yaxis=dict(showline=True, linewidth=1, linecolor='black', mirror=True),
            plot_bgcolor='rgba(0,0,0,0)'
        )

        return fig
    
    @app.callback(
        Output('mapping-pie', 'figure'),
        Input('interval', 'n_intervals'),
        Input('mapping-checklist', 'value')
    )
    def update_pie(_, selected_categories):
        """
        Update the pie chart

        """
        if not os.path.exists(flagstats_csv):
            return go.Figure()  # returing blank figure
        
        # Load data
        df = pd.read_csv(flagstats_csv)

        # Pie values
        pie_values = df[selected_categories].sum().tolist()

        # Plot data
        plot_data = [
            go.Pie(labels=selected_categories,
                   values=pie_values,
                   marker=dict(colors=[COL_MAP[cat] for cat in selected_categories]))
        ]
        fig = go.Figure(data=plot_data)

        return fig


    @app.callback(
        Output('mapping-barplot', 'figure'),
        Input('interval', 'n_intervals'),
        Input('mapping-checklist', 'value')
    )
    def update_barplot(_, selected_categories):
        """
        Update the mapping barplot

        """

        if not os.path.exists(flagstats_csv):
            return go.Figure()  # returing blank figure
        
        # Load data
        df = pd.read_csv(flagstats_csv)

        # Prepare plotting information
        plot_data = [
            go.Bar(x=df["barcode"],
                   y=df[cat],
                   marker=dict(color=COL_MAP[cat]),
                   name=cat)
            for cat in selected_categories
        ]

        # Generate figure
        fig = go.Figure(data=plot_data)

        # Format
        fig.update_layout(
            #title="Read Mapping Statistics",
            yaxis_title="Num. of Reads",

            barmode="stack",
            
            xaxis=dict(showline=True, linewidth=1, linecolor='black', mirror=True),
            yaxis=dict(showline=True, linewidth=1, linecolor='black', mirror=True),
            
            plot_bgcolor="rgba(0,0,0,0)",
            hovermode="x"
        )
        fig.update_traces(
            marker=dict(line=dict(color="black", width=1))
        )

        return fig

    @app.callback(
        Output('expt-summary', 'children'),
        Input('interval', 'n_intervals')
    )
    def update_expt_summary(_):
        """
        Update the experimental summary text

        t0.strftime("%Y-%m-%d %H:%M:%S")

        """

        t1 = datetime.datetime.now().replace(microsecond=0)
        text = f"Experiment: {expt_name}<br>"
        text += f"Started at: {t0.strftime('%Y-%m-%d %H:%M:%S')}<br>"
        children = [
            html.H3("Run Details"),
            html.P([
            f"Experiment: {expt_name}",
            html.Br(),
            f"Started at: {t0.strftime('%Y-%m-%d %H:%M:%S')}",
            html.Br(),
            f"Time elapsed: {t1 - t0}"])
        ]
        return children

    app.layout = html.Div(id="overall", 
                          children=layout)

    return app


# For testing purposes
if __name__ == "__main__":
    
    expt_name = "0000-00-00_example"
    flagstats_csv = f"results/{expt_name}/summary.bam_flagstats.csv"
    bedcov_csv = f"results/{expt_name}/summary.bedcov.csv"

    dashboard = create_dashboard(
        expt_name,
        flagstats_csv,
        bedcov_csv
    )

    dashboard.run_server(debug=True)
