import logging
import threading
from abc import ABC, abstractmethod
from dash import Dash, html, dcc

from nomadic.util.regions import RegionBEDParser
from nomadic.realtime.dashboard.components import (
    ExperimentSummaryFASTQ,
    OverallGauge,
    MappingStatsPie,
    MappingStatsBarplot,
    RegionCoveragePie,
    RegionCoverageStrip,
    DepthProfileLinePlot,
    DepthProfileHistogram,
)


# --------------------------------------------------------------------------------
# Abstract base class
#
# --------------------------------------------------------------------------------


class RealtimeDashboardBuilder(ABC):

    FAST_NAME = "fast-interval"
    FAST_INTERVAL = 1_000
    SLOW_NAME = "interval"
    SLOW_INTERVAL = 60_000

    def __init__(self, expt_name, css_style_sheets):
        self.expt_name = expt_name
        self.css_style_sheets = css_style_sheets

        self.components = []

    def _gen_app(self):
        """
        Generate an instance of a Dash app

        """

        app = Dash(name=__name__, external_stylesheets=self.css_style_sheets)
        app_log = logging.getLogger("werkzeug")
        app_log.setLevel(logging.ERROR)

        return app

    @abstractmethod
    def _gen_layout(self):
        """
        Define the layout of the dashboard

        This will link with the stylesheet that is used to produce
        the overall dashboard organisation

        """
        pass

    def _gen_timer(self, name, speed):
        """
        Generate the timer which is a feature of all real-time dashboards

        """

        return dcc.Interval(id=name, interval=speed, n_intervals=0)

    def run(self, in_thread: bool = False, **kwargs):
        """
        Run the dashboard

        """

        app = self._gen_app()
        layout = self._gen_layout()
        layout.append(self._gen_timer(self.FAST_NAME, self.FAST_INTERVAL))
        layout.append(self._gen_timer(self.SLOW_NAME, self.SLOW_INTERVAL))

        for component in self.components:
            component.callback(app)

        app.layout = html.Div(id="overall", children=layout)

        if in_thread:
            dashboard_thread = threading.Thread(target=app.run, name="dashboard")
            dashboard_thread.start()
        else:
            app.run(**kwargs)


# --------------------------------------------------------------------------------
# Concrete dashboards
#
# --------------------------------------------------------------------------------


class MappingRTDashboard(RealtimeDashboardBuilder):
    """
    Build a dashboard with a focus on mapping statistics

    """

    # TODO: this is bad, need better way to handle repeated constants
    CATEGORIES = ["n_uniq_mapped", "n_chim_mapped", "n_mult_mapped", "n_unmapped"]
    CSS_STYLES = ["assets/mapping-style.css"]

    def __init__(
        self,
        expt_name,
        regions: RegionBEDParser,
        fastq_csv,
        flagstats_csv,
        bedcov_csv,
        depth_csv,
    ):
        """
        Initialise all of the dashboard components

        """

        super().__init__(expt_name, self.CSS_STYLES)

        self.regions = regions
        self.fastq_csv = fastq_csv
        self.flagstats_csv = flagstats_csv
        self.bedcov_csv = bedcov_csv
        self.depth_csv = depth_csv

        # Initialise all the components
        self.expt_summary = ExperimentSummaryFASTQ(
            expt_name=expt_name, component_id="expt-summary", fastq_csv=fastq_csv
        )
        self.mapping_pie = MappingStatsPie(
            expt_name,
            component_id="mapping-pie",
            flagstats_csv=flagstats_csv,
            checklist_id="mapping-checklist",
        )
        self.mapping_bar = MappingStatsBarplot(
            expt_name,
            component_id="mapping-barplot",
            flagstats_csv=flagstats_csv,
            checklist_id="mapping-checklist",
        )

        self.region_pie = RegionCoveragePie(
            expt_name,
            component_id="bedcov-pie",
            regions=regions,
            bedcov_csv=bedcov_csv,
            dropdown_id="bedcov-dropdown",
        )
        self.region_strip = RegionCoverageStrip(
            expt_name,
            component_id="bedcov-strip",
            regions=regions,
            bedcov_csv=bedcov_csv,
            dropdown_id="bedcov-dropdown",
        )
        self.depth_line = DepthProfileLinePlot(
            expt_name,
            regions=self.regions,
            component_id="depth-line",
            depth_csv=self.depth_csv,
            region_dropdown_id="depth-dropdown",
        )
        self.depth_hist = DepthProfileHistogram(
            expt_name,
            regions=self.regions,
            component_id="depth-hist",
            depth_csv=self.depth_csv,
            region_dropdown_id="depth-dropdown",
        )

        # Put them into the components
        self.components.append(self.expt_summary)
        self.components.append(self.mapping_pie)
        self.components.append(self.mapping_bar)
        self.components.append(self.region_pie)
        self.components.append(self.region_strip)
        self.components.append(self.depth_line)
        self.components.append(self.depth_hist)

    def _gen_layout(self):
        """
        Generate the layout

        """

        # Input elements
        checklist = dcc.Checklist(
            id="mapping-checklist",
            options=self.CATEGORIES,
            value=self.CATEGORIES,
            inline=True,
        )

        dropdown = dcc.Dropdown(
            id="bedcov-dropdown",
            options=["mean_cov", "n_reads", "cov_gr100_per"],
            value="n_reads",
            style=dict(width="300px"),
        )

        depth_dropdown = dcc.Dropdown(
            id="depth-dropdown",
            options=self.regions.names,
            value="kelch13",
            style=dict(width="300px"),
        )

        # Layout
        banner = html.Div(className="banner", children=self.expt_summary.get_layout())

        mapping_row = html.Div(
            className="mapping-row",
            children=[
                html.H3("Read Mapping Statistics", style=dict(marginTop="0px")),
                checklist,
                html.Div(
                    className="mapping-plots",
                    children=[
                        self.mapping_pie.get_layout(),
                        self.mapping_bar.get_layout(),
                    ],
                ),
            ],
        )

        bedcov_row = html.Div(
            className="bedcov-row",
            children=[
                html.H3("Region Coverage Statistics", style=dict(marginTop="0px")),
                dropdown,
                html.Div(
                    className="bedcov-plots",
                    children=[
                        self.region_pie.get_layout(),
                        self.region_strip.get_layout(),
                    ],
                ),
            ],
        )

        depth_row = html.Div(
            className="depth-row",
            children=[
                html.H3("Region Coverage Profiles", style=dict(marginTop="0px")),
                depth_dropdown,
                html.Div(
                    className="depth-plots",
                    children=[
                        self.depth_hist.get_layout(),
                        self.depth_line.get_layout(),
                    ],
                ),
            ],
        )

        layout = [banner, mapping_row, bedcov_row, depth_row]

        return layout
