import logging
import threading
from abc import ABC, abstractmethod
from dash import Dash, html, dcc

from nomadic.util.regions import RegionBEDParser
from nomadic.util.metadata import MetadataTableParser
from nomadic.realtime.dashboard.components import (
    ExperimentSummaryFASTQ,
    MappingStatsPie,
    MappingStatsBarplot,
    RegionCoveragePie,
    RegionCoverageStrip,
    DepthProfileLinePlot,
    DepthProfileHistogram,
    DepthProfileCumulativeDist,
    VariantHeatmap
)


# --------------------------------------------------------------------------------
# Interface for different dashboards
#
# TODO:
# - There is a lot of duplication in subclasses; I can probably replicate
#   what was done in the experiment pipeline heirarchy, placing reused elements
#   in the superclass (e.g. def _add_banner(), def _add_depth_row()) to reduce this.
# - Let's do this while I am here...
# --------------------------------------------------------------------------------



class RealtimeDashboardBuilder(ABC):
    """
    Interface for all real-time dashboards

    TODO: Very confused about whether I need to `self.` things
    for them to persist

    """

    FAST_NAME = "fast-interval"
    FAST_INTERVAL = 1_000
    SLOW_NAME = "interval"
    SLOW_INTERVAL = 60_000

    def __init__(self, expt_name, css_style_sheets):
        self.expt_name = expt_name
        self.css_style_sheets = css_style_sheets

        self.components = []
        self.layout = []

    @abstractmethod
    def _gen_layout(self):
        """
        Define the layout of the dashboard

        This will link with the stylesheet that is used to produce
        the overall dashboard organisation

        """
        pass

    def _gen_app(self):
        """
        Generate an instance of a Dash app

        """

        app = Dash(name=__name__, external_stylesheets=self.css_style_sheets)
        app_log = logging.getLogger("werkzeug")
        app_log.setLevel(logging.ERROR)

        return app

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
        self._gen_layout()
        self.layout.append(self._gen_timer(self.FAST_NAME, self.FAST_INTERVAL))
        self.layout.append(self._gen_timer(self.SLOW_NAME, self.SLOW_INTERVAL))

        for component in self.components:
            component.callback(app)

        app.layout = html.Div(id="overall", children=self.layout)

        if in_thread:
            dashboard_thread = threading.Thread(target=app.run, name="dashboard")
            dashboard_thread.start()
        else:
            app.run(**kwargs)

    # ---------------------------------------------------------------------------
    # Below here, we are putting concrete methods for creating different
    # pieces of a dashboard that may be shared across dashboard subclasse
    #
    # ---------------------------------------------------------------------------

    def _add_banner(self, fastq_csv: str) -> None:
        """
        Add a banner which shows the logo and summarise the number of FASTQ files
        processed

        """
        
        # Create the component
        self.expt_summary = ExperimentSummaryFASTQ(
            expt_name=self.expt_name, 
            component_id="expt-summary", 
            fastq_csv=fastq_csv
        )

        # Define banner layout
        banner = html.Div(className="banner", children=self.expt_summary.get_layout())
        
        # Add to components and layout
        self.components.append(self.expt_summary)
        self.layout.append(banner)


    def _add_mapping_row(self, flagstats_csv: str) -> None:
        """
        Add a row summarising mapping statistics for all barcodes

        """
        
        # Here you can set what mapping statistics to make viewable
        CATEGORIES = ["n_primary", "n_chimeria", "n_secondary", "n_unmapped"]

        # Create checklist of mapping categories
        checklist = dcc.Checklist(
            id="mapping-checklist",
            options=CATEGORIES,
            value=CATEGORIES,
            inline=True,
        )

        # Create a pieplot component
        self.mapping_pie = MappingStatsPie(
            self.expt_name,
            component_id="mapping-pie",
            flagstats_csv=flagstats_csv,
            checklist_id="mapping-checklist",
        )

        # Create a barplot component
        self.mapping_bar = MappingStatsBarplot(
            self.expt_name,
            component_id="mapping-barplot",
            flagstats_csv=flagstats_csv,
            checklist_id="mapping-checklist",
        )

        # Define the row layout
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

        # Add to components and layout
        self.components.append(self.mapping_pie)
        self.components.append(self.mapping_bar)
        self.layout.append(mapping_row)

    def _add_bedcov_row(self, bedcov_csv: str, regions: RegionBEDParser) -> None:
        """
        Add a row summarising coverage across amplicons
        
        """

        # Dropdown of different statistics
        dropdown = dcc.Dropdown(
            id="bedcov-dropdown",
            options=["mean_cov", "n_reads", "cov_gr100_per"],
            value="n_reads",
            style=dict(width="300px"),
        )

        # Define components
        self.region_pie = RegionCoveragePie(
            self.expt_name,
            component_id="bedcov-pie",
            regions=regions,
            bedcov_csv=bedcov_csv,
            dropdown_id="bedcov-dropdown",
        )
        self.region_strip = RegionCoverageStrip(
            self.expt_name,
            component_id="bedcov-strip",
            regions=regions,
            bedcov_csv=bedcov_csv,
            dropdown_id="bedcov-dropdown",
        )

        # Define the row layout
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

        # Add to components and layout
        self.components.append(self.region_pie)
        self.components.append(self.region_strip)
        self.layout.append(bedcov_row)

    def _add_depth_row(self, depth_csv: str, regions: RegionBEDParser) -> None:
        """
        Add a row summarising coverage across amplicons
        
        """

        depth_dropdown = dcc.Dropdown(
            id="depth-dropdown",
            options=regions.names,
            value=regions.names[0],
            #value="kelch13",
            style=dict(width="300px"),
        )

        self.depth_line = DepthProfileLinePlot(
            self.expt_name,
            regions=regions,
            component_id="depth-line",
            depth_csv=depth_csv,
            region_dropdown_id="depth-dropdown",
        )
        self.depth_hist = DepthProfileCumulativeDist(
            self.expt_name,
            regions=regions,
            component_id="depth-hist",
            depth_csv=depth_csv,
            region_dropdown_id="depth-dropdown",
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

        # Add to components and layout
        self.components.append(self.depth_line)
        self.components.append(self.depth_hist)
        self.layout.append(depth_row)


# --------------------------------------------------------------------------------
# Concrete dashboards
#
# --------------------------------------------------------------------------------


class MappingRTDashboard(RealtimeDashboardBuilder):
    """
    Build a dashboard with a focus on mapping statistics

    """

    CSS_STYLE = ["assets/mapping-style.css"]

    def __init__(
        self,
        expt_name,
        regions: RegionBEDParser,
        metadata: MetadataTableParser,
        fastq_csv,
        flagstats_csv,
        bedcov_csv,
        depth_csv
    ):
        """
        Initialise all of the dashboard components

        """

        super().__init__(expt_name, self.CSS_STYLE)

        self.regions = regions
        self.metadata = metadata
        self.fastq_csv = fastq_csv
        self.flagstats_csv = flagstats_csv
        self.bedcov_csv = bedcov_csv
        self.depth_csv = depth_csv

    def _gen_layout(self):
        """
        Generate the layout

        """
        self._add_banner(self.fastq_csv)
        self._add_mapping_row(self.flagstats_csv)
        self._add_bedcov_row(self.bedcov_csv, self.regions)
        self._add_depth_row(self.depth_csv, self.regions)


class CallingRTDashboard(RealtimeDashboardBuilder):
    """
    Build a dashboard with a focus on mapping statistics and
    variant calling results

    """

    CSS_STYLE = ["assets/calling-style.css"]

    def __init__(
        self,
        expt_name,
        regions: RegionBEDParser,
        metadata: MetadataTableParser,
        fastq_csv,
        flagstats_csv,
        bedcov_csv,
        depth_csv,
        variant_csv
    ):
        """
        Initialise all of the dashboard components

        """

        super().__init__(expt_name, self.CSS_STYLE)

        self.regions = regions
        self.metadata = metadata
        self.fastq_csv = fastq_csv
        self.flagstats_csv = flagstats_csv
        self.bedcov_csv = bedcov_csv
        self.depth_csv = depth_csv
        self.variant_csv = variant_csv

    def _add_calling_row(self, variant_csv: str, regions: RegionBEDParser) -> None:
        """
        Add a row summarising variant calling results across amplicons
        
        """

        # Dropdown of different target regions
        variant_dropdown = dcc.Dropdown(
            id="variant-dropdown",
            options=self.regions.names,
            value=regions.names[0],
            style=dict(width="300px"),
        )

        # Define component
        self.variant_heat = VariantHeatmap(
            self.expt_name,
            regions=regions,
            metadata=self.metadata,
            component_id="variant-heat",
            variant_csv=self.variant_csv,
            region_dropdown_id="variant-dropdown",
        )

        # Define layout
        variant_row = html.Div(
            className="variant-row",
            children=[
                html.H3("Preliminary Variant Calling", style=dict(marginTop="0px")),
                variant_dropdown,
                html.Div(
                    className="variant-plots",
                    children=[
                        self.variant_heat.get_layout()
                    ],
                ),
            ],
        )

        # Add to components and layout
        self.components.append(self.variant_heat)
        self.layout.append(variant_row)

    def _gen_layout(self):
        """
        Generate the layout for the variant calling dashboard

        """
        self._add_banner(self.fastq_csv)
        self._add_mapping_row(self.flagstats_csv)
        self._add_bedcov_row(self.bedcov_csv, self.regions)
        self._add_depth_row(self.depth_csv, self.regions)
        self._add_calling_row(self.variant_csv, self.regions)

