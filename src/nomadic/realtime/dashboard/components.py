import os
import datetime
import pandas as pd

from abc import ABC, abstractmethod
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output

import numpy as np
import seaborn as sns
import plotly.graph_objects as go
from matplotlib.colors import rgb2hex

from nomadic.util.regions import RegionBEDParser


# --------------------------------------------------------------------------------
# Parameters
#
# --------------------------------------------------------------------------------


TIMER_INTERVAL_ID = "interval"

MAPPING_CATS = ["n_uniq_mapped", "n_chim_mapped", "n_mult_mapped", "n_unmapped"]
MAPPING_COLS = dict(
    zip(
        MAPPING_CATS,
        [
            rgb2hex(c)
            for c in sns.color_palette("Blues_r", len(MAPPING_CATS) - 1)
            + [(0.75, 0.75, 0.75)]
        ],
    )
)


# --------------------------------------------------------------------------------
# Interface for a single real-time dashboard component
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

# --------------------------------------------------------------------------------
# OVERALL STATISTICS
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
            Output(self.component_id, "children"), Input("interval", "n_intervals")
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


class ExperimentSummaryFASTQ(RealtimeDashboardComponent):
    """
    Overview of the experimental state, and number of FASTQ files
    processed

    """

    logo_src_path = "assets/nomadic_logo-01.png"

    def __init__(self, expt_name: str, component_id: str, fastq_csv: str):
        super().__init__(expt_name, component_id)
        self.t0 = datetime.datetime.now().replace(microsecond=0)
        self.fastq_csv = fastq_csv

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
            Output(self.component_id, "children"), Input("fast-interval", "n_intervals")
        )
        def _update(_):
            """Called every time an input changes"""

            # Re-load the number of FASTQ processed
            n_fastq = 0
            if os.path.exists(self.fastq_csv):
                df = pd.read_csv(self.fastq_csv)
                n_fastq = df["n_processed_fastq"].sum()

            # Update time
            t1 = datetime.datetime.now().replace(microsecond=0)

            # Define the text section
            # TODO: this is pretty horrible to make look decent
            tab = "\t"
            n_tabs = 4
            children = [
                html.H3("Run Overview", style=dict(margin="0px", marginTop="20px")),
                html.Pre(
                    # pre_contents,
                    [
                        f"Experiment Name:{tab*(n_tabs-1)}{self.expt_name}",
                        html.Br(),
                        f"Started at:{tab*n_tabs}{self.t0.strftime('%Y-%m-%d %H:%M:%S')}",
                        html.Br(),
                        f"Time elapsed:{tab*n_tabs}{t1 - self.t0}",
                        html.Br(),
                        f"No. FASTQ Processed:{tab*(n_tabs-2)}{n_fastq}",
                    ],
                    style=dict(fontFamily="Arial", margin="0px"),
                ),
            ]

            return children


# --------------------------------------------------------------------------------
# MAPPING STATISTICS
# --------------------------------------------------------------------------------


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

            MAR = 20
            fig.update_layout(margin=dict(t=MAR, l=MAR, r=MAR, b=MAR), showlegend=False)

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
                yaxis=dict(
                    showline=True,
                    linewidth=1,
                    linecolor="black",
                    mirror=True,
                    showgrid=True,
                    gridcolor="lightgray",
                    gridwidth=0.5,
                    griddash="dot",
                ),
                plot_bgcolor="rgba(0,0,0,0)",
                hovermode="x",
                legend=dict(
                    orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0
                ),
            )
            fig.update_traces(marker=dict(line=dict(color="black", width=1)))

            return fig


# --------------------------------------------------------------------------------
# BED COVERAGE SUMMARIES
# --------------------------------------------------------------------------------


class RegionCoveragePie(RealtimeDashboardComponent):
    """
    Make a pie chart that shows read mapping statistics

    """

    def __init__(
        self,
        expt_name: str,
        component_id: str,
        regions: RegionBEDParser,
        bedcov_csv: str,
        dropdown_id: str,
    ):

        # Store inputs
        super().__init__(expt_name, component_id)
        self.bedcov_csv = bedcov_csv
        self.regions = regions
        self.dropdown_id = dropdown_id

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
            Input(self.dropdown_id, "value"),
        )
        def _update(_, dropdown_stat):
            """Called every time an input changes"""

            # Load data
            if not os.path.exists(self.bedcov_csv):
                return go.Figure()
            df = pd.read_csv(self.bedcov_csv)
            df["name"] = pd.Categorical(
                values=df["name"], categories=self.regions.names, ordered=True
            )

            # Compute totals
            pie_data = df.groupby("name")[dropdown_stat].sum()

            # Generate figure
            fig = go.Figure(
                data=[
                    go.Pie(
                        values=pie_data.values,
                        labels=pie_data.index,
                        marker=dict(
                            colors=[
                                self.regions.col_map_hex[name]
                                for name in pie_data.index
                            ]
                        ),
                        sort=False,
                    )
                ]
            )

            MAR = 20
            fig.update_layout(showlegend=False, margin=dict(t=MAR, l=MAR, r=MAR, b=MAR))

            return fig


class RegionCoverageStrip(RealtimeDashboardComponent):
    """
    Make a stripplot that shows read mapping statistics

    """

    def __init__(
        self,
        expt_name: str,
        component_id: str,
        regions: RegionBEDParser,
        bedcov_csv: str,
        dropdown_id: str,
    ):

        # Store inputs
        super().__init__(expt_name, component_id)
        self.bedcov_csv = bedcov_csv
        self.regions = regions
        self.dropdown_id = dropdown_id

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
            Input(self.dropdown_id, "value"),
        )
        def _update(_, dropdown_stat, shift=0):
            """Called every time an input changes"""

            # Load data, and sort
            if not os.path.exists(self.bedcov_csv):
                return go.Figure()
            df = pd.read_csv(self.bedcov_csv)
            df["name"] = pd.Categorical(
                values=df["name"], categories=self.regions.names, ordered=True
            )

            # Prepare plotting data
            plot_data = [
                go.Scatter(
                    x=tdf[
                        "barcode"
                    ],  # + np.random.uniform(-shift, shift, tdf.shape[0]),
                    y=tdf[dropdown_stat],
                    mode="markers",
                    marker=dict(size=10, color=self.regions.col_map_hex[target]),
                    name=target,
                )
                for target, tdf in df.groupby("name")
            ]

            # Fix y-axis minimum at zero
            min_y = 0
            max_y = df[dropdown_stat].max() * 1.1

            # Create plot
            fig = go.Figure()
            for plot_trace in plot_data:
                fig.add_trace(plot_trace)

            fig.update_layout(
                yaxis_title=dropdown_stat,
                xaxis=dict(showline=True, linewidth=1, linecolor="black", mirror=True),
                yaxis=dict(
                    showline=True,
                    linewidth=1,
                    linecolor="black",
                    mirror=True,
                    showgrid=True,
                    gridcolor="lightgray",
                    gridwidth=0.5,
                    griddash="dot",
                    range=[min_y, max_y],
                ),
                plot_bgcolor="rgba(0,0,0,0)",
                legend=dict(
                    orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0
                ),
            )

            return fig


class OverallGauge(RealtimeDashboardComponent):
    """
    Present a gauge indicating what percentage of all
    basepairs have exceeded the set coverage threshold

    """

    def __init__(self, expt_name: str, component_id: str, bedcov_csv: str):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.bedcov_csv = bedcov_csv

    def _define_layout(self):
        """
        The layout for this graph
        """
        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        """
        Define the update callback for the gauge
        """

        @app.callback(
            Output(self.component_id, "figure"), Input(self.interval_id, "n_intervals")
        )
        def _update(_):
            """Update each interval"""

            # Load data, and sort
            if not os.path.exists(self.bedcov_csv):
                return go.Figure()
            df = pd.read_csv(self.bedcov_csv)

            # Compute key statistics
            total_bp = df["length"].sum()
            covered_bp = df["cov_gr100"].sum()
            per_covered = 100 * (covered_bp / total_bp)

            # Create the figure
            fig = go.Figure(
                go.Indicator(
                    domain={"x": [0, 1], "y": [0, 1]},
                    value=per_covered,
                    # gauge = {'shape': "bullet"},
                    number=dict(suffix="%"),
                    mode="gauge+number+delta",
                    title={"text": ""},
                    # delta = {'reference': 94, 'suffix': '%'},
                    gauge={
                        "axis": {"range": [None, 100]},
                        "steps": [
                            {"range": [0, 50], "color": "lightgrey"},
                        ],
                        "threshold": {
                            "line": {"color": "steelblue", "width": 4},
                            "thickness": 1,
                            "value": 95,
                        },
                    },
                    #'shape': 'bullet'}
                )
            )

            # Add text annotation
            MAR = 10
            fig.update_layout(
                font=dict(size=8),
                margin=dict(t=MAR, l=MAR, r=MAR, b=MAR),
                annotations=[
                    go.layout.Annotation(
                        x=0.5,
                        y=-0.15,
                        xref="paper",
                        yref="paper",
                        text="<B>Target Coverage Achieved (of all bp)</B>",
                        showarrow=False,
                        font=dict(size=10, color="black", family="Helvetica"),
                    )
                ],
            )

            return fig


# --------------------------------------------------------------------------------
# DEPTH SUMMARIES
# --------------------------------------------------------------------------------


class DepthProfileLinePlot(RealtimeDashboardComponent):
    """
    Make a depth profile plot

    """

    pal = "inferno"

    def __init__(
        self,
        expt_name: str,
        regions: RegionBEDParser,
        component_id: str,
        depth_csv: str,
        region_dropdown_id: str,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.regions = regions
        self.depth_csv = depth_csv
        self.region_dropdown_id = region_dropdown_id

    def _define_layout(self):
        """Layout is graph"""
        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        """Check the timer, selected region, and plot all lines"""

        @app.callback(
            Output(self.component_id, "figure"),
            Input(self.interval_id, "n_intervals"),
            Input(self.region_dropdown_id, "value"),
        )
        def _update(_, target_region):
            """Called every time an input changes"""

            # Load data
            if not os.path.exists(self.depth_csv):
                return go.Figure()
            df = pd.read_csv(self.depth_csv)

            # Select target gene and group by barcode
            region_df = df.query("name == @target_region")
            grps = region_df.groupby("barcode")
            n_barcodes = len(grps)
            col_map = dict(
                zip(
                    grps.groups.keys(),
                    [rgb2hex(c) for c in sns.color_palette(self.pal, n_barcodes)],
                )
            )

            # Control axes
            min_y = 0
            max_y = region_df["depth"].max() * 1.1

            # Plot data
            plot_data = [
                go.Scatter(
                    x=bdf["pos"],
                    y=bdf["depth"],
                    mode="lines",
                    marker=dict(color=col_map[barcode_name]),
                    name=barcode_name,
                )
                for barcode_name, bdf in region_df.groupby("barcode")
            ]

            # Create the plot
            fig = go.Figure(plot_data)

            MAR = 40
            # Format
            fig.update_layout(
                title=dict(text=f"{target_region}"),
                margin=dict(t=MAR, l=MAR, r=MAR, b=MAR),
                yaxis_title="Depth",
                xaxis_title="Genomic Position",
                xaxis=dict(showline=True, linewidth=1, linecolor="black", mirror=True),
                yaxis=dict(
                    showline=True,
                    linewidth=1,
                    linecolor="black",
                    mirror=True,
                    showgrid=True,
                    gridcolor="lightgray",
                    gridwidth=0.5,
                    griddash="dot",
                    range=[min_y, max_y],
                ),
                plot_bgcolor="rgba(0,0,0,0)",
                hovermode="x",
                legend=dict(
                    orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=0
                ),
                showlegend=False,
            )

            return fig


class DepthProfileHistogram(RealtimeDashboardComponent):
    """
    Make a histogram of the depth profile

    """

    pal = "inferno"

    def __init__(
        self,
        expt_name: str,
        regions: RegionBEDParser,
        component_id: str,
        depth_csv: str,
        region_dropdown_id: str,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.regions = regions
        self.depth_csv = depth_csv
        self.region_dropdown_id = region_dropdown_id

    def _define_layout(self):
        """Layout is graph"""
        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        """Check the timer, selected region, and plot all lines"""

        @app.callback(
            Output(self.component_id, "figure"),
            Input(self.interval_id, "n_intervals"),
            Input(self.region_dropdown_id, "value"),
        )
        def _update(_, target_region):
            """Called every time an input changes"""

            # Load data
            if not os.path.exists(self.depth_csv):
                return go.Figure()
            df = pd.read_csv(self.depth_csv)

            # Select target gene and group by barcode
            region_df = df.query("name == @target_region")
            grps = region_df.groupby("barcode")
            n_barcodes = len(grps)
            col_map = dict(
                zip(
                    grps.groups.keys(),
                    [rgb2hex(c) for c in sns.color_palette(self.pal, n_barcodes)],
                )
            )

            # Control axes
            min_x = 0
            max_x = region_df["depth"].max() * 1.1

            # Plot data
            plot_data = [
                go.Histogram(x=bdf["depth"], marker=dict(color=col_map[barcode_name]), name=barcode_name)
                for barcode_name, bdf in region_df.groupby("barcode")
            ]

            # Create the plot
            fig = go.Figure(plot_data)

            MAR = 40
            # Format
            fig.update_layout(
                barmode="stack",
                margin=dict(t=MAR, l=MAR, r=MAR, b=MAR),
                yaxis_title="Count",
                xaxis_title="Depth",
                xaxis=dict(
                    showline=True,
                    linewidth=1,
                    linecolor="black",
                    mirror=True,
                    range=[min_x, max_x],
                ),
                yaxis=dict(
                    showline=True,
                    linewidth=1,
                    linecolor="black",
                    mirror=True,
                    showgrid=True,
                    gridcolor="lightgray",
                    gridwidth=0.5,
                    griddash="dot",
                ),
                plot_bgcolor="rgba(0,0,0,0)",
                hovermode="x",
                legend=dict(
                    orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=0
                ),
                showlegend=False,
            )

            return fig
