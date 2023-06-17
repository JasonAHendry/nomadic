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
)


# --------------------------------------------------------------------------------
# Abstract base class
#
# --------------------------------------------------------------------------------


class RealtimeDashboardBuilder(ABC):
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

    def _gen_timer(self):
        """
        Generate the timer which is a feature of all real-time dashboards

        """

        return dcc.Interval(id="interval", interval=1_000, n_intervals=0)

    def run(self, in_thread: bool = False, **kwargs):
        """
        Run the dashboard

        """

        app = self._gen_app()
        layout = self._gen_layout()
        layout.append(self._gen_timer())

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


class TestingDashboard(RealtimeDashboardBuilder):
    """
    Build a basic dashboard with this styling

    """

    CATEGORIES = ["n_uniq_mapped", "n_chim_mapped", "n_mult_mapped", "n_unmapped"]

    def __init__(self, expt_name, css_style_sheets, flagstats_csv):
        """Testing this architecture"""
        super().__init__(expt_name, css_style_sheets)

        self.flagstats_csv = flagstats_csv

        # Initialise all the components
        self.mapping_pie = MappingStatsPie(
            expt_name,
            component_id="mapping-pie",
            flagstats_csv=flagstats_csv,
            checklist_id="mapping-checklist",
        )
        self.components.append(self.mapping_pie)

    def _gen_layout(self):
        """
        Generate the layout

        """
        checklist = dcc.Checklist(
            id="mapping-checklist",
            options=self.CATEGORIES,
            value=self.CATEGORIES,
            inline=True,
        )

        layout = [
            html.H2("Testing"),
            html.Hr(),
            self.mapping_pie.get_layout(),
            checklist,
        ]

        return layout


class MappingRTDashboard(RealtimeDashboardBuilder):
    """
    Build a dashboard with a focus on mapping statistics

    """

    # TODO: this is bad, need better way to handle repeated constants
    CATEGORIES = ["n_uniq_mapped", "n_chim_mapped", "n_mult_mapped", "n_unmapped"]

    def __init__(
        self,
        expt_name,
        regions: RegionBEDParser,
        css_style_sheets,
        fastq_csv,
        flagstats_csv,
        bedcov_csv,
    ):
        """
        Initialise all of the dashboard components

        """

        super().__init__(expt_name, css_style_sheets)

        self.regions = regions
        self.fastq_csv = fastq_csv
        self.flagstats_csv = flagstats_csv
        self.bedcov_csv = bedcov_csv

        # Initialise all the components
        self.expt_summary = ExperimentSummaryFASTQ(
            expt_name=expt_name, component_id="expt-summary", fastq_csv=fastq_csv
        )
        # self.gauge = OverallGauge(expt_name,
        #     "gauge",
        #     bedcov_csv
        # )

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

        # Put them into the components
        self.components.append(self.expt_summary)
        #self.components.append(self.gauge)
        self.components.append(self.mapping_pie)
        self.components.append(self.mapping_bar)
        self.components.append(self.region_pie)
        self.components.append(self.region_strip)

    def _gen_layout(self):
        """
        Generate the layout

        """

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
            style=dict(width="300px")
        )

        top_row = html.Div(className="top-row", 
                           children=self.expt_summary.get_layout())
        
        middle_row = html.Div(
            className="middle-row",
            children=[                
                html.H3("Read Mapping Statistics", style=dict(marginTop="0px")),
                checklist,
                html.Div(className='mapping-plots',
                         children=[self.mapping_pie.get_layout(), 
                                   self.mapping_bar.get_layout()]),
                ]
        )

        bottom_row = html.Div(
            className="bottom-row",
            children=[
                html.H3("Region Coverage Statistics", style=dict(marginTop="0px")),
                dropdown,
                html.Div(className='region-plots',  
                         children=[self.region_pie.get_layout(), 
                                   self.region_strip.get_layout()])
                ]
        )

        layout = [
            top_row,
            middle_row,
            bottom_row,
            html.Br(),
            html.Br(),
        ]

        return layout
