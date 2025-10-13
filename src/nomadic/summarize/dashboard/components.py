# import datetime
# import os
from abc import ABC, abstractmethod
# from typing import Optional
# import re

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

# --------------------------------------------------------------------------------
# Interface for a single real-time dashboard component
#
# --------------------------------------------------------------------------------


class SummaryDashboardComponent(ABC):
    """
    Interface for summary dashboard components

    """

    def __init__(self, summary_name: str, component_id: str):
        """
        Initialise components of the dashboard component

        """

        # User
        self.summary_name = summary_name
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


class ThroughputSummary(SummaryDashboardComponent):
    """
    Make a pie chart that shows read mapping statistics

    """

    logo_src_path = "assets/nomadic_logo.png"

    def __init__(self, summary_name: str, throughput_csv: str, component_id: str):
        self.throughput_csv = throughput_csv
        self.throughput_df = pd.read_csv(throughput_csv, index_col="sample_type")
        super().__init__(summary_name, component_id)

    def _define_layout(self):
        """
        Define the layout to be a dcc.Graph object with the
        appropriate ID

        """

        layout = html.Div(
            className="logo-and-summary",
            children=[
                html.Img(id="logo", src=self.logo_src_path),
                html.Div(
                    id="throughput-summary",
                    children=[
                        html.H3("Throughput Details"),
                        html.P(
                            [
                                f"Experiments: {self.throughput_df.columns.shape[0] - 1}",
                                html.Br(),
                                f"Field samples (total): {self.throughput_df.loc['field', 'All']}",
                                html.Br(),
                                f"Field samples (unique): {self.throughput_df.loc['field_unique', 'All']}",
                                html.Br(),
                            ]
                        ),
                    ],
                ),
            ],
        )

        return layout

    def callback(self, app: Dash) -> None:
        """
        Define the update callback for the pie chart

        """


class SamplesPie(SummaryDashboardComponent):
    """
    Make a pie chart that shows read mapping statistics

    """

    def __init__(
        self,
        summary_name: str,
        samples_csv: str,
        component_id: str,
    ):
        self.samples_csv = samples_csv
        self.df = pd.read_csv(samples_csv)
        self.df = self.df.groupby("status").count()["sample_id"]
        super().__init__(summary_name, component_id)

    def _define_layout(self):
        """
        Define the layout to be a dcc.Graph object with the
        appropriate ID

        """
        fig = go.Figure(
            data=[
                go.Pie(
                    values=self.df.values,
                    labels=self.df.index,
                    sort=False,
                    hole=0.3,
                )
            ]
        )

        MAR = 20
        fig.update_layout(showlegend=False, margin=dict(t=MAR, l=MAR, r=MAR, b=MAR))

        return dcc.Graph(id=self.component_id, figure=fig)

    def callback(self, app: Dash) -> None:
        """
        Define the update callback for the pie chart

        """


class QualityControl(SummaryDashboardComponent):
    STATISTICS = [
        "mean_cov_field",
        "per_field_passing",
        "per_field_contam",
        "per_field_lowcov",
    ]

    def __init__(
        self, summary_name: str, coverage_csv: str, component_id: str, dropdown_id: str
    ) -> None:
        """
        Initialisation loads the coverage data and prepares for plotting;

        """

        self.coverage_csv = coverage_csv
        self.coverage_df = pd.read_csv(coverage_csv)
        self.plot_df = pd.pivot_table(
            index="expt_name",
            columns="name",
            values=self.STATISTICS,
            dropna=False,
            observed=False,
            data=self.coverage_df,
        )

        self.dropdown_id = dropdown_id
        super().__init__(summary_name, component_id)

    def _define_layout(self):
        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        @app.callback(
            Output(self.component_id, "figure"), Input(self.dropdown_id, "value")
        )
        def _update(focus_stat: str):
            """Called whenver the input changes"""
            plot_data = [
                go.Heatmap(
                    x=self.plot_df[focus_stat].columns,
                    y=self.plot_df[focus_stat].index,
                    z=self.plot_df[focus_stat],
                    text=self.plot_df[focus_stat],
                    # colorscale="Reds",
                    xgap=1,
                    ygap=1,
                    colorbar=dict(
                        title=focus_stat, outlinecolor="black", outlinewidth=1
                    ),
                    hoverongaps=False,
                    # **STAT_KWARGS[STAT]
                )
            ]
            MAR = 40
            fig = go.Figure(plot_data)
            fig.update_layout(
                width=1200,
                height=600,
                yaxis_title="Experiments",
                hovermode="y unified",
                paper_bgcolor="white",  # Sets the background color of the paper
                plot_bgcolor="white",
                title=dict(text=focus_stat),
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
            fig.update_traces(
                text=self.plot_df[focus_stat],
                texttemplate="%{text:.0f}",
                textfont_size=12,
            )
            return fig


class PrevalenceBarplot(SummaryDashboardComponent):
    GENE_SETS = {
        "Resistance": ["crt", "dhps", "dhfr", "kelch13", "mdr1"],
        "Diversity": ["ama1", "csp"],
    }

    def __init__(
        self,
        summary_name: str,
        prevalence_csv: str,
        component_id: str,
        radio_id: str,
    ) -> None:
        """
        Initialisation loads the coverage data and prepares for plotting;

        """

        self.prevalence_csv = prevalence_csv
        self.prev_df = pd.read_csv(prevalence_csv)

        self.radio_id = radio_id
        super().__init__(summary_name, component_id)

    def _define_layout(self):
        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        @app.callback(
            Output(self.component_id, "figure"), Input(self.radio_id, "value")
        )
        def _update(gene_set: str):
            """Called whenver the input changes"""

            genes = self.GENE_SETS[gene_set]

            # Limit to key genes
            plot_df = self.prev_df.query("gene in @genes")
            plot_df.sort_values(["gene", "chrom", "pos"], inplace=True)

            # Prepare plotting data
            customdata = np.stack(
                [
                    plot_df["n_samples"],
                    plot_df["n_passed"],
                    plot_df["n_mixed"] + plot_df["n_mut"],
                ],
                axis=-1,
            )

            # Plotting
            htemp = "%{y:0.1f}% (%{customdata[2]}/%{customdata[1]})"
            plot_data = [
                go.Bar(
                    x=plot_df["mutation"],
                    y=plot_df["prevalence"],
                    customdata=customdata,
                    hovertemplate=htemp,
                    name="Prevalence",
                    error_y=dict(
                        type="data",
                        array=plot_df["prevalence_highci"] - plot_df["prevalence"],
                        arrayminus=plot_df["prevalence"] - plot_df["prevalence_lowci"],
                    ),
                ),
                # go.Bar(
                #     x=plot_df["mutation"],
                #     y=plot_df["per_mixed"],
                #     customdata=customdata,
                #     hovertemplate=htemp,
                #     name="Mixed",
                # ),
            ]
            fig = go.Figure(plot_data)
            fig.update_layout(
                yaxis_title="Prevalence (%)",
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
                barmode="stack",
                legend=dict(
                    orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0
                ),
                plot_bgcolor="rgba(0,0,0,0)",
                hovermode="x unified",
            )
            fig.update_yaxes(range=[0, 100])
            fig.update_traces(marker=dict(line=dict(color="black", width=1)))

            return fig
