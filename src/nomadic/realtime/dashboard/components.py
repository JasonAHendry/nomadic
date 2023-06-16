import os
import datetime
import pandas as pd

from abc import ABC, abstractmethod
from dash import Dash, html, dcc
from dash.dependencies import Input, Output

import seaborn as sns
import plotly.graph_objects as go
from matplotlib.colors import rgb2hex


# --------------------------------------------------------------------------------
# Parameters
#
# --------------------------------------------------------------------------------


TIMER_INTERVAL_ID = "interval"

MAPPING_CATS = ["n_uniq_mapped", "n_chim_mapped", "n_mult_mapped", "n_unmapped"]
MAPPING_COLS = dict(
    zip(
        MAPPING_CATS,
        [rgb2hex(c) for c in sns.color_palette("Blues_r") + [(0.75, 0.75, 0.75)]],
    )
)


# --------------------------------------------------------------------------------
# Abstract base class
#
# --------------------------------------------------------------------------------


class RealtimeDashboardComponent(ABC):
    """
    Interface for real-time dashboard components

    These respond to a timer, which is defined
    by the `interval_id`.

    """

    interval_id = TIMER_INTERVAL_ID

    def __init__(self, expt_name: str, component_id: str):
        """
        Initialise components of the dashboard component

        """

        # User
        self.expt_name = expt_name
        self.component_id = component_id

        # Always
        self.layout_obj = self._define_layout()

    @abstractmethod
    def _define_layout(self):
        pass

    def get_layout(self):
        """
        Get the object that will be passed to the overall HTML
        layout

        Typically, this is something from `html.` or `dcc.`, e.g.
        `dcc.Graph`

        """

        if self.layout_obj is None:
            raise ValueError("Must define `self.layout_component`.")

        return self.layout_obj

    @abstractmethod
    def callback(self, app: Dash) -> None:
        """
        Define the callback for this componenet, which will cause it to update
        in response to the timer, as well (potentially) other inputs

        """
        pass


# --------------------------------------------------------------------------------
# Concrete components
#
# --------------------------------------------------------------------------------


class ExperimentSummary(RealtimeDashboardComponent):
    """
    Make a pie chart that shows read mapping statistics

    """

    logo_src_path = "assets/nomadic_logo-01.png"

    def __init__(self, expt_name: str, component_id: str):
        super().__init__(expt_name, component_id)
        self.t0 = datetime.datetime.now().replace(microsecond=0)

    def _define_layout(self):
        """
        Define the layout to be a dcc.Graph object with the
        appropriate ID

        """

        layout = html.Div(
            className="logo-and-summary",
            children=[
                html.Img(id="logo", src=self.logo_src_path),
                html.Div(id="expt-summary"),
            ],
        )

        return layout

    def callback(self, app: Dash) -> None:
        """
        Define the update callback for the pie chart

        """

        @app.callback(
            Output(self.component_id, "children"), 
            Input("interval", "n_intervals")
        )
        def _update(_):
            """Called every time an input changes"""

            t1 = datetime.datetime.now().replace(microsecond=0)

            children = [
                html.H3("Run Details"),
                html.P(
                    [
                        f"Experiment: {self.expt_name}",
                        html.Br(),
                        f"Started at: {self.t0.strftime('%Y-%m-%d %H:%M:%S')}",
                        html.Br(),
                        f"Time elapsed: {t1 - self.t0}",
                    ]
                ),
            ]

            return children


class MappingStatsPie(RealtimeDashboardComponent):
    """
    Make a pie chart that shows read mapping statistics

    """

    def __init__(
        self, expt_name: str, component_id: str, flagstats_csv: str, checklist_id: str
    ):

        # Store inputs
        super().__init__(expt_name, component_id)
        self.flagstats_csv = flagstats_csv
        self.checklist_id = checklist_id

    def _define_layout(self):
        """
        Define the layout to be a dcc.Graph object with the
        appropriate ID

        """

        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        """
        Define the update callback for the pie chart

        """

        @app.callback(
            Output(self.component_id, "figure"),
            Input(self.interval_id, "n_intervals"),
            Input(self.checklist_id, "value"),
        )
        def _update(_, selected_categories):
            """Called every time an input changes"""

            # Load data
            if not os.path.exists(self.flagstats_csv):
                return go.Figure()
            df = pd.read_csv(self.flagstats_csv)

            # Compute totals
            pie_values = df[selected_categories].sum().tolist()

            # Generate figure
            fig = go.Figure(
                data=[
                    go.Pie(
                        labels=selected_categories,
                        values=pie_values,
                        marker=dict(
                            colors=[MAPPING_COLS[cat] for cat in selected_categories]
                        ),
                    )
                ]
            )

            return fig


class MappingStatsBarplot(RealtimeDashboardComponent):
    """
    Make a bar chart that shows read mapping statistics

    """

    def __init__(
        self, expt_name: str, component_id: str, flagstats_csv: str, checklist_id: str
    ):

        # Store inputs
        super().__init__(expt_name, component_id)
        self.flagstats_csv = flagstats_csv
        self.checklist_id = checklist_id

    def _define_layout(self):
        """
        Define the layout to be a dcc.Graph object with the
        appropriate ID

        """

        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        """
        Define the update callback for the bar chart

        """

        @app.callback(
            Output(self.component_id, "figure"),
            Input(self.interval_id, "n_intervals"),
            Input(self.checklist_id, "value"),
        )
        def _update(_, selected_categories):
            """Called every time an input changes"""

            # Load data
            if not os.path.exists(self.flagstats_csv):
                return go.Figure()
            df = pd.read_csv(self.flagstats_csv)

            # Generate figure
            fig = go.Figure(
                data=[
                    go.Bar(
                        x=df["barcode"],
                        y=df[cat],
                        marker=dict(color=MAPPING_COLS[cat]),
                        name=cat,
                    )
                    for cat in selected_categories
                ]
            )

            # Format
            fig.update_layout(
                yaxis_title="No. Reads",
                barmode="stack",
                xaxis=dict(showline=True, linewidth=1, linecolor="black", mirror=True),
                yaxis=dict(showline=True, linewidth=1, linecolor="black", mirror=True),
                plot_bgcolor="rgba(0,0,0,0)",
                hovermode="x",
            )
            fig.update_traces(marker=dict(line=dict(color="black", width=1)))

            return fig
