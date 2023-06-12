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
    app = Dash(__name__)
    app_log = logging.getLogger("werkzeug")
    app_log.setLevel(logging.ERROR)

    layout = [
        html.Img(src="assets/nomadic_logo-01.png", style={"width": "400px"}),
        html.Div(id="expt-summary", style={"margin": "20px", "fontFamily": "Helvetica"}),
        dcc.Graph(id='mapping-barplot'),
        dcc.Checklist(id='mapping-checklist',
                      options=CATEGORIES,
                      value=CATEGORIES,
                      inline=True),
        dcc.Interval(
            id='interval',
            interval=1_000,
            n_intervals=0
        )
    ]

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
            title="Read Mapping Statistics",
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
        element = html.P([
            f"Experiment: {expt_name}",
            html.Br(),
            f"Started at: {t0.strftime('%Y-%m-%d %H:%M:%S')}",
            html.Br(),
            f"Time elapsed: {t1 - t0}"
        ])
        return element

    app.layout = html.Div(layout, style={"fontFamily": "Helvetica", "margin": "20px"})

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
