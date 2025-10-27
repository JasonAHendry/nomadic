import datetime
import os
from abc import ABC, abstractmethod
from typing import Optional
import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from matplotlib.colors import rgb2hex

from nomadic.util.metadata import MetadataTableParser
from nomadic.util.regions import RegionBEDParser
from i18n import t

pd.options.mode.chained_assignment = None


# --------------------------------------------------------------------------------
# Parameters
#
# --------------------------------------------------------------------------------


TIMER_INTERVAL_ID = "interval"

MAPPING_CATS = ["n_primary", "n_supplementary", "n_secondary", "n_unmapped"]
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

    logo_src_path = "assets/nomadic_logo.png"

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

    logo_src_path = "assets/nomadic_logo.png"

    def __init__(
        self,
        expt_name: str,
        component_id: str,
        fastq_csv: str,
        start_time: Optional[datetime.datetime] = None,
        is_realtime: bool = True,
    ):
        super().__init__(expt_name, component_id)
        if start_time is not None:
            self.t0 = start_time
        else:
            self.t0 = datetime.datetime.now().replace(microsecond=0)
        self.fastq_csv = fastq_csv
        self.is_realtime = is_realtime

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

            content = [
                f"{t('Experiment Name')}:{tab * (n_tabs - 1)}{self.expt_name}",
                html.Br(),
                f"{t('Started at')}:{tab * n_tabs}{self.t0.strftime('%Y-%m-%d %H:%M:%S')}",
                html.Br(),
                f"{t('Time elapsed')}:{tab * n_tabs}{t1 - self.t0}",
                html.Br(),
                f"{t('FASTQs Processed')}:{tab * (n_tabs - 2)}{n_fastq}",
            ]

            if not self.is_realtime:
                # remove time elapsed
                del content[-3]
                del content[-2]

            children = [
                html.H3("Run Overview", style=dict(margin="0px", marginTop="20px")),
                html.Pre(
                    # pre_contents,
                    content,
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
        self,
        expt_name: str,
        component_id: str,
        read_mapping_csv: str,
        checklist_id: str,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.read_mapping_csv = read_mapping_csv
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
            if not os.path.exists(self.read_mapping_csv):
                return go.Figure()
            df = pd.read_csv(self.read_mapping_csv)

            if "n_chimeria" in df.columns:
                # for backwards compatibility to be able to show old experiments where this old name was used
                df.rename(columns={"n_chimeria": "n_supplementary"}, inplace=True)

            # Compute totals
            pie_cats = [c for c in MAPPING_CATS if c in selected_categories]
            pie_values = df[pie_cats].sum().tolist()

            # Generate figure
            fig = go.Figure(
                data=[
                    go.Pie(
                        labels=[t(category) for category in selected_categories],
                        values=pie_values,
                        marker=dict(colors=[MAPPING_COLS[cat] for cat in pie_cats]),
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
        self,
        expt_name: str,
        component_id: str,
        read_mapping_csv: str,
        checklist_id: str,
        metadata: MetadataTableParser,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.read_mapping_csv = read_mapping_csv
        self.checklist_id = checklist_id
        self.metadata = metadata

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
            if not os.path.exists(self.read_mapping_csv):
                return go.Figure()
            df = pd.read_csv(self.read_mapping_csv)

            if "n_chimeria" in df.columns:
                # for backwards compatibility to be able to show old experiments where this old name was used
                df.rename(columns={"n_chimeria": "n_supplementary"}, inplace=True)

            if "sample_id" not in df.columns:
                # for backwards compatibility to be able to show old experiments where this column was not in the data
                df = df.join(self.metadata.sample_ids_df, on="barcode")

            x = df.apply(
                sample_string_from_row,
                axis=1,
            )

            # Generate figure
            fig = go.Figure(
                data=[
                    go.Bar(
                        x=x,
                        y=df[cat],
                        marker=dict(color=MAPPING_COLS[cat]),
                        name=t(cat),
                    )
                    for cat in MAPPING_CATS
                    if cat in selected_categories
                ]
            )

            # Format
            fig.update_layout(
                xaxis_title="Samples",
                yaxis_title="Alignment Count",
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
                hovermode="x unified",
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
        region_coverage_csv: str,
        dropdown_id: str,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.region_coverage_csv = region_coverage_csv
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
            if not os.path.exists(self.region_coverage_csv):
                return go.Figure()
            df = pd.read_csv(self.region_coverage_csv)
            df["name"] = pd.Categorical(
                values=df["name"], categories=self.regions.names, ordered=True
            )

            # Compute totals
            pie_data = df.groupby("name", observed=False)[dropdown_stat].sum()

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
        region_coverage_csv: str,
        dropdown_id: str,
        metadata: MetadataTableParser,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.region_coverage_csv = region_coverage_csv
        self.regions = regions
        self.dropdown_id = dropdown_id
        self.metadata = metadata

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
            if not os.path.exists(self.region_coverage_csv):
                return go.Figure()
            df = pd.read_csv(self.region_coverage_csv)
            df["name"] = pd.Categorical(
                values=df["name"], categories=self.regions.names, ordered=True
            )
            if "sample_id" not in df.columns:
                # for backwards compatibility to be able to show old experiments where this column was not in the data
                df = df.join(self.metadata.sample_ids_df, on="barcode")

            # Prepare plotting data
            # plot_data = [
            #     go.Scatter(
            #         x=tdf[
            #             "barcode"
            #         ],  # + np.random.uniform(-shift, shift, tdf.shape[0]),
            #         y=tdf[dropdown_stat],
            #         mode="markers",
            #         marker=dict(size=10, color=self.regions.col_map_hex[target]),
            #         name=target,
            #     )
            #     for target, tdf in df.groupby("name")
            # ]

            # Fix y-axis minimum at zero
            min_y = 0
            max_y = df[dropdown_stat].max() * 1.1

            # Create plot
            # fig = go.Figure()
            # for plot_trace in plot_data:
            #     fig.add_trace(plot_trace)
            plot_df = pd.DataFrame.from_dict(
                {
                    t("Sample"): df.apply(sample_string_from_row, axis=1),
                    t("name"): df["name"],
                    t(dropdown_stat): df[dropdown_stat],
                }
            )
            fig = px.strip(
                plot_df,
                x=t("Sample"),
                color=t("name"),
                color_discrete_map=self.regions.col_map_hex,
                y=t(dropdown_stat),
            )

            fig.update_layout(
                xaxis_title=t("Samples"),
                yaxis_title=t(dropdown_stat),
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
            fig.update_traces(marker=dict(size=10))

            return fig


class OverallGauge(RealtimeDashboardComponent):
    """
    Present a gauge indicating what percentage of all
    basepairs have exceeded the set coverage threshold

    """

    def __init__(self, expt_name: str, component_id: str, region_coverage_csv: str):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.region_coverage_csv = region_coverage_csv

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
            if not os.path.exists(self.region_coverage_csv):
                return go.Figure()
            df = pd.read_csv(self.region_coverage_csv)

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

    pal = "tab20"

    def __init__(
        self,
        expt_name: str,
        regions: RegionBEDParser,
        component_id: str,
        depth_profiles_csv: str,
        region_dropdown_id: str,
        metadata: MetadataTableParser,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.regions = regions
        self.depth_profiles_csv = depth_profiles_csv
        self.region_dropdown_id = region_dropdown_id
        self.metadata = metadata

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
            if not os.path.exists(self.depth_profiles_csv):
                return go.Figure()
            df = pd.read_csv(self.depth_profiles_csv)

            # Select target gene and group by barcode
            region_df = df.query("name == @target_region")

            grps = region_df.groupby("barcode", observed=False)
            n_samples = len(grps)
            col_map = dict(
                zip(
                    grps.groups.keys(),
                    [rgb2hex(c) for c in sns.color_palette(self.pal, n_samples)],
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
                    marker=dict(color=col_map[barcode]),
                    name=sample_string(barcode, self.metadata.get_sample_id(barcode)),
                )
                for barcode, bdf in grps
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
                hovermode="x unified",
                legend=dict(
                    orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=0
                ),
                showlegend=False,
            )

            return fig


class DepthProfileCumulativeDist(RealtimeDashboardComponent):
    """
    Make a cumulative distribution plot of depth across
    a given region

    """

    pal = "tab20"

    def __init__(
        self,
        expt_name: str,
        regions: RegionBEDParser,
        component_id: str,
        depth_profiles_csv: str,
        region_dropdown_id: str,
        metadata: MetadataTableParser,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.regions = regions
        self.depth_profiles_csv = depth_profiles_csv
        self.region_dropdown_id = region_dropdown_id
        self.metadata = metadata

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
            if not os.path.exists(self.depth_profiles_csv):
                return go.Figure()
            df = pd.read_csv(self.depth_profiles_csv)

            # Select target gene and group by barcode
            region_df = df.query("name == @target_region")
            grps = region_df.groupby("barcode", observed=False)
            n_samples = len(grps)
            col_map = dict(
                zip(
                    grps.groups.keys(),
                    [rgb2hex(c) for c in sns.color_palette(self.pal, n_samples)],
                )
            )

            target_info = self.regions.df.query("name == @target_region").squeeze()
            target_bp = target_info["end"] - target_info["start"]

            # Plot data
            percentiles = 100 * (1 - np.arange(target_bp + 1) / target_bp)
            plot_data = [
                go.Scatter(
                    x=sorted(sample_df["depth"]),
                    y=percentiles,
                    marker=dict(color=col_map[barcode]),
                    mode="lines",
                    name=sample_string(barcode, self.metadata.get_sample_id(barcode)),
                )
                for barcode, sample_df in grps
            ]

            # Create the plot
            fig = go.Figure(plot_data)

            MAR = 40
            # Format
            fig.update_layout(
                barmode="stack",
                margin=dict(t=MAR, l=MAR, r=MAR, b=MAR),
                yaxis_title="Region Covered (%)",
                xaxis_title="Minimum Depth",
                xaxis=dict(
                    showline=True,
                    linewidth=1,
                    linecolor="black",
                    mirror=True,
                    # range=[min_x, max_x],
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
                hovermode="x unified",
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
        depth_profiles_csv: str,
        region_dropdown_id: str,
        metadata: MetadataTableParser,
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.regions = regions
        self.depth_profiles_csv = depth_profiles_csv
        self.region_dropdown_id = region_dropdown_id
        self.metadata = metadata

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
            if not os.path.exists(self.depth_profiles_csv):
                return go.Figure()
            df = pd.read_csv(self.depth_profiles_csv)

            # Select target gene and group by sample
            region_df = df.query("name == @target_region")
            grps = region_df.groupby("barcode", observed=False)
            n_samples = len(grps)
            col_map = dict(
                zip(
                    grps.groups.keys(),
                    [rgb2hex(c) for c in sns.color_palette(self.pal, n_samples)],
                )
            )

            # Control axes
            min_x = 0
            max_x = region_df["depth"].max() * 1.1

            # Plot data
            plot_data = [
                go.Histogram(
                    x=bdf["depth"],
                    marker=dict(color=col_map[barcode]),
                    name=sample_string(barcode, self.metadata.get_sample_id(barcode)),
                )
                for barcode, bdf in grps
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


# --------------------------------------------------------------------------------
# VARIANT CALLING SUMMARIES
# --------------------------------------------------------------------------------


class VariantHeatmap(RealtimeDashboardComponent):
    """
    Make a heatmap of variant calls

    """

    MUT_SET = ["synonymous", "missense"]

    def __init__(
        self,
        expt_name: str,
        metadata: MetadataTableParser,
        regions: RegionBEDParser,
        component_id: str,
        variant_csv: str,
        region_dropdown_id: str,
        known_variant_csv: str = "",  # should be optional
    ):
        # Store inputs
        super().__init__(expt_name, component_id)
        self.regions = regions
        self.metadata = metadata
        self.component_id = component_id
        barcodes = (
            metadata.barcodes.copy()
        )  # TODO: would this be by reference or value?
        if "unclassified" in barcodes:
            barcodes.remove("unclassified")
        self.categories = [
            sample_string(barcode, self.metadata.get_sample_id(barcode))
            for barcode in barcodes
        ]

        self.variant_csv = variant_csv
        self.region_dropdown_id = region_dropdown_id

        if known_variant_csv:
            self.known = True
            self.known_df = pd.read_csv(known_variant_csv)

    def _define_layout(self):
        """Layout is graph"""
        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        @app.callback(
            Output(self.component_id, "figure"),
            Input(self.interval_id, "n_intervals"),
            Input(self.region_dropdown_id, "value"),
        )
        def _update(_, target_region):
            """Called every time an input changes"""

            # Load data
            if not os.path.exists(self.variant_csv):
                return go.Figure()
            df = pd.read_csv(self.variant_csv)

            # Filter on target, variant types, and passing QC
            qry = (
                "amplicon == @target_region"
                " and mut_type in @self.MUT_SET"
                " and gt != './.'"
            )
            target_df = df.query(qry)
            if "sample_id" not in target_df.columns:
                # for backwards compatibility to be able to show old experiments where this column was not in the data
                target_df = target_df.join(self.metadata.sample_ids_df, on="barcode")

            # Munge for plot
            target_df["sample_string"] = pd.Categorical(
                values=target_df.apply(sample_string_from_row, axis=1),
                categories=self.categories,
                ordered=True,
            )

            # Pivot
            plot_df = pd.pivot_table(
                index="aa_change",
                columns="sample_string",
                values=["wsaf", "dp", "gt"],
                aggfunc=lambda x: x,
                data=target_df,
                dropna=False,
                observed=False,
            )

            # Reorder mutations based on position
            mutations = target_df[["pos", "aa_change"]].drop_duplicates("aa_change")
            mutations.sort_values(["pos", "aa_change"], inplace=True)
            plot_df.index = pd.Categorical(
                values=plot_df.index, categories=mutations.aa_change, ordered=True
            )
            plot_df.sort_index(inplace=True)

            # Hover statment
            customdata = np.stack([plot_df["dp"], plot_df["gt"]], axis=-1)
            htemp = "<b>%{x}</b><br>"
            htemp += "<b>WSAF:</b> %{z:0.3f}<br>"
            htemp += "<b>Depth:</b> %{customdata[0]}<br>"
            htemp += "<b>Genotype:</b> %{customdata[1]}<br>"

            # Plot
            plot_data = [
                go.Heatmap(
                    x=plot_df["wsaf"].columns,
                    y=plot_df["wsaf"].index,
                    z=plot_df["wsaf"],
                    customdata=customdata,
                    zmin=0,
                    zmax=1,
                    xgap=1,
                    ygap=1,
                    colorscale="Spectral_r",
                    colorbar=dict(title="WSAF*", outlinecolor="black", outlinewidth=1),
                    hovertemplate=htemp,
                    hoverongaps=False,
                    name="",
                )
            ]

            fig = go.Figure(plot_data)

            # Format
            MAR = 40
            # SZ = 50
            fig.update_layout(
                xaxis_title="Samples",
                hovermode="y unified",
                paper_bgcolor="white",  # Sets the background color of the paper
                plot_bgcolor="white",
                title=dict(text=target_region),
                margin=dict(t=MAR, l=MAR, r=MAR, b=MAR),
                xaxis=dict(
                    showline=True, linecolor="black", linewidth=2, dtick=1, mirror=True
                ),
                yaxis=dict(
                    showline=True, linecolor="black", linewidth=2, dtick=1, mirror=True
                ),
                xaxis_showgrid=False,
                yaxis_showgrid=False,
                # height=n_mutations*SZ # TOOD: how to adjust dynamically
            )

            return fig


def sample_string_from_row(row: pd.Series) -> str:
    """Create the representation of sample"""

    if "sample_id" not in row or not isinstance(row["sample_id"], str):
        return sample_string(row["barcode"])

    return sample_string(row["barcode"], row["sample_id"])


def sample_string(barcode: str, sample_id: Optional[str] = None) -> str:
    if sample_id is None:
        return barcode

    if barcode == "unclassified":
        return barcode

    match = re.search(r"\d+", barcode)

    assert match, "should have found a number"

    number = int(match.group())

    return f"{sample_id} ({number:02d})"
