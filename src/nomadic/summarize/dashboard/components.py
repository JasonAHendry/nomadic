from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, dcc, html
from dash.dependencies import Input, Output

from nomadic.summarize.compute import (
    compute_variant_prevalence,
    compute_variant_prevalence_per,
)
from i18n import t

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


SAMPLE_COLORS = {
    "missing": "#636EFA",
    "failing": "#EF553B",
    "passing": "#00CC96",
}


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
        self.n = len(self.df)
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
                    textinfo="label+percent+value",
                    marker=dict(colors=[SAMPLE_COLORS[cat] for cat in self.df.index]),
                )
            ]
        )

        MAR = 20
        fig.update_layout(
            showlegend=True,
            margin=dict(t=MAR, l=MAR, r=MAR, b=MAR),
            annotations=[
                dict(
                    text=f"N={self.n}", font_size=20, showarrow=False, xanchor="center"
                )
            ],
        )

        return dcc.Graph(id=self.component_id, figure=fig)

    def callback(self, app: Dash) -> None:
        """
        Define the update callback for the pie chart

        """


COV_MAX = 10000
COLORSCALES = {
    "per_field_passing": [
        [0.00, "#FF0033"],
        [0.50, "#FF9900"],
        [0.70, "#FFEB33"],  # set 70% as okay
        [0.90, "#80FF7E"],  # above 90% is good
        [1.00, "#3A9A3E"],
    ],
    "per_field_contam": [
        [0.00, "#3A9A3E"],  # low contamination is good
        [0.20, "#FFEB33"],  # set 30% as okay
        [0.30, "#FF9900"],  # contamination is bad
        [1.00, "#FF0033"],
    ],
    "per_field_lowcov": [
        [0.00, "#3A9A3E"],  # low lowcov is good
        [0.10, "#80FF7E"],  # above 90% is good
        [0.30, "#FFEB33"],  # set 30% as okay
        [0.50, "#FF9900"],  # lowcov is bad
        [1.00, "#FF0033"],
    ],
    "mean_cov_field": [
        [0, "#FF0033"],  # low coverage is bad
        [25 / COV_MAX, "#FF9900"],
        [50 / COV_MAX, "#FFEB33"],  # set threshold as okay
        [200 / COV_MAX, "#80FF7E"],  # above 200 we can make good calls
        [500 / COV_MAX, "#3A9A3E"],  # above 500 is excellent also for low freq calls
        [1.0, "#7585FE"],  # coverage is uncapped
    ],
}


class AmpliconsBarplot(SummaryDashboardComponent):
    """
    Make a bar chart that shows the Amplicons Statistics

    """

    def __init__(
        self,
        summary_name: str,
        component_id: str,
        samples_amplicons_csv: str,
    ):
        df = pd.read_csv(samples_amplicons_csv)
        # Store inputs
        plot_df = pd.crosstab(df["name"], df["status"])
        n_missing = (df["status"] == "missing").sum()
        missing = [n_missing] * len(plot_df.index)
        # Generate figure
        fig = go.Figure(
            data=[
                go.Bar(
                    x=plot_df.index,
                    y=missing if column == "missing" else plot_df[column],
                    texttemplate="%{y}",
                    name=column,
                    marker=dict(color=SAMPLE_COLORS[column]),
                )
                for column in ["passing", "failing", "missing"]
            ]
        )

        # Format
        fig.update_layout(
            xaxis_title="Amplicons",
            yaxis_title="Number of Samples",
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
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        )
        fig.update_traces(marker=dict(line=dict(color="black", width=1)))
        self.fig = fig
        super().__init__(summary_name, component_id)

    def _define_layout(self):
        """
        Define the layout to be a dcc.Graph object with the
        appropriate ID

        """

        return dcc.Graph(id=self.component_id, figure=self.fig)

    def callback(self, app: Dash) -> None:
        pass


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
            legend = (
                "Field samples (%)"
                if "per_" in focus_stat
                else "Mean Coverage"
                if "cov" in focus_stat
                else ""
            )
            plot_data = [
                go.Heatmap(
                    x=self.plot_df[focus_stat].columns,
                    y=self.plot_df[focus_stat].index,
                    z=self.plot_df[focus_stat],
                    text=self.plot_df[focus_stat],
                    xgap=1,
                    ygap=1,
                    zmin=0,
                    zmax=100 if "per_" in focus_stat else 10000,
                    colorbar=dict(title=legend, outlinecolor="black", outlinewidth=1),
                    hoverongaps=False,
                    colorscale=COLORSCALES[focus_stat],
                )
            ]
            MAR = 40
            fig = go.Figure(plot_data)
            fig.update_layout(
                width=1200,
                height=600,
                xaxis_title="Amplicons",
                yaxis_title="Experiments",
                hovermode="y unified",
                paper_bgcolor="white",  # Sets the background color of the paper
                plot_bgcolor="white",
                title=dict(text=t(focus_stat)),
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
            unit = "%" if "per_" in focus_stat else "x" if "cov" in focus_stat else ""
            fig.update_traces(
                text=self.plot_df[focus_stat],
                texttemplate="%{text:.0f}" + unit,
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
        analysis_csv: str,
        master_csv: str,
        component_id: str,
        radio_id: str,
        radio_id_by: str,
    ) -> None:
        """
        Initialisation loads the coverage data and prepares for plotting;

        """

        self.analysis_csv = analysis_csv
        self.analysis_df = pd.read_csv(analysis_csv)

        self.master__csv = master_csv
        self.master_df = pd.read_csv(master_csv)

        self.radio_id = radio_id
        self.radio_id_by = radio_id_by
        super().__init__(summary_name, component_id)

    def _define_layout(self):
        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        @app.callback(
            Output(self.component_id, "figure"),
            Input(self.radio_id, "value"),
            Input(self.radio_id_by, "value"),
        )
        def _update(gene_set: str, by: str):
            """Called whenver the input changes"""

            genes = self.GENE_SETS[gene_set]  # noqa: F841 later used in query

            # Limit to key genes
            analysis_df = self.analysis_df.query("gene in @genes")
            if by == "All":
                plot_df = compute_variant_prevalence(analysis_df)
            else:
                plot_df = compute_variant_prevalence_per(
                    analysis_df, self.master_df, by.split("_")
                )
                if "_" in by:
                    # we need to create this column
                    plot_df[by] = (
                        plot_df[by.split("_")].astype(str).agg("_".join, axis=1)
                    )
            plot_df.sort_values(["gene", "chrom", "pos"], inplace=True)

            data = []
            htemp = "%{y:0.1f}% (%{customdata[2]}/%{customdata[1]})"

            if by == "All":
                # Prepare plotting data
                customdata = np.stack(
                    [
                        plot_df["n_samples"],
                        plot_df["n_passed"],
                        plot_df["n_mixed"] + plot_df["n_mut"],
                    ],
                    axis=-1,
                )
                data.append(
                    go.Bar(
                        x=plot_df["mutation"],
                        y=plot_df["prevalence"],
                        customdata=customdata,
                        hovertemplate=htemp,
                        name="Prevalence",
                        error_y=dict(
                            type="data",
                            array=plot_df["prevalence_highci"] - plot_df["prevalence"],
                            arrayminus=plot_df["prevalence"]
                            - plot_df["prevalence_lowci"],
                        ),
                    )
                )
            else:
                for group in plot_df[by].unique():
                    group_df = plot_df.query(f"{by} == @group")
                    # Prepare plotting data
                    customdata = np.stack(
                        [
                            group_df["n_samples"],
                            group_df["n_passed"],
                            group_df["n_mixed"] + group_df["n_mut"],
                        ],
                        axis=-1,
                    )
                    data.append(
                        go.Bar(
                            x=group_df["mutation"],
                            y=group_df["prevalence"],
                            customdata=customdata,
                            hovertemplate=htemp,
                            name=str(group),
                            error_y=dict(
                                type="data",
                                array=plot_df["prevalence_highci"]
                                - plot_df["prevalence"],
                                arrayminus=plot_df["prevalence"]
                                - plot_df["prevalence_lowci"],
                            ),
                        )
                    )

            # Plotting
            fig = go.Figure(data)
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
                legend=dict(
                    orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0
                ),
                plot_bgcolor="rgba(0,0,0,0)",
                hovermode="x unified",
            )
            fig.update_yaxes(range=[0, 100])
            fig.update_traces(marker=dict(line=dict(color="black", width=1)))

            return fig


class PrevalenceHeatmap(SummaryDashboardComponent):
    """
    Make a heatmap of prevalences

    """

    def __init__(
        self,
        summary_name: str,
        prevalence_region_csv: str,
        component_id: str,
        gene_dropdown_id: str,
    ):
        self.gene_dropdown_id = gene_dropdown_id
        self.df = pd.read_csv(prevalence_region_csv)
        super().__init__(summary_name, component_id)

    def _define_layout(self):
        """Layout is graph"""
        return dcc.Graph(id=self.component_id)

    def callback(self, app: Dash) -> None:
        @app.callback(
            Output(self.component_id, "figure"),
            Input(self.gene_dropdown_id, "value"),
        )
        def _update(target_gene):
            """Called every time an input changes"""

            df = self.df.query("gene == @target_gene")
            plot_df = pd.pivot_table(
                index="aa_change",
                columns="region",
                values=["prevalence", "n_mixed", "n_mut", "n_passed"],
                data=df,
            )

            # Sort, so aa changes are in correct order
            pos_order = (
                df.drop_duplicates("aa_change")
                .set_index("aa_change")["aa_pos"]
                .sort_values(ascending=True)
                .index
            )
            plot_df = plot_df.reindex(pos_order)

            # Hover statment
            customdata = np.stack(
                [plot_df["n_mixed"], plot_df["n_mut"], plot_df["n_passed"]], axis=-1
            )
            htemp = "<b>%{y} (%{x})</b><br>"
            htemp += "<b>Prevalence:</b> %{z:.0f}%<br>"
            htemp += "<b>Samples:</b> %{customdata[2]}<br>"
            htemp += "<b>Mixed:</b> %{customdata[0]}<br>"
            htemp += "<b>Clonal:</b> %{customdata[1]}<br>"

            plot_data = [
                go.Heatmap(
                    x=plot_df["prevalence"].columns,
                    y=plot_df["prevalence"].index,
                    z=plot_df["prevalence"],
                    texttemplate="%{z:.0f}%",
                    customdata=customdata,
                    zmin=0,
                    zmax=100,
                    xgap=1,
                    ygap=1,
                    colorscale="Spectral_r",
                    colorbar=dict(title="", outlinecolor="black", outlinewidth=1),
                    hoverongaps=False,
                    hovertemplate=htemp,
                    name="",
                )
            ]
            MAR = 40
            fig = go.Figure(plot_data)
            fig.update_layout(
                xaxis_title="Regions",
                hovermode="y unified",
                paper_bgcolor="white",  # Sets the background color of the paper
                plot_bgcolor="white",
                title=dict(text=target_gene),
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
