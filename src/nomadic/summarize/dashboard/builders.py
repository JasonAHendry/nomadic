import logging
import threading
from abc import ABC, abstractmethod
from dash import Dash, html, dcc
from datetime import datetime
from typing import Optional

from i18n import t

# from importlib.resources import files, as_file
from nomadic.summarize.dashboard.components import (
    AmpliconsBarplot,
    PrevalenceHeatmap,
    SamplesPie,
    ThroughputSummary,
    QualityControl,
    PrevalenceBarplot,
)


class SummaryDashboardBuilder(ABC):
    """
    Interface for summary dashboards

    """

    def __init__(self, summary_name, css_style_sheets):
        self.summary_name = summary_name
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

    # def _gen_timer(self, name, speed):
    #     """
    #     Generate the timer which is a feature of all real-time dashboards

    #     """

    #     return dcc.Interval(id=name, interval=speed, n_intervals=0)

    def run(self, in_thread: bool = False, **kwargs):
        """
        Run the dashboard

        """

        # setup_translations()
        app = self._gen_app()
        self._gen_layout()

        for component in self.components:
            component.callback(app)

        app.layout = html.Div(id="overall", children=self.layout)

        if in_thread:
            dashboard_thread = threading.Thread(
                target=lambda: app.run(**kwargs), name="dashboard", daemon=True
            )
            dashboard_thread.start()
        else:
            app.run(**kwargs)

    # ---------------------------------------------------------------------------
    # Below here, we are putting concrete methods for creating different
    # pieces of a dashboard that may be shared across dashboard subclasse
    #
    # ---------------------------------------------------------------------------

    def _add_throughput_banner(self, throughput_csv: str) -> None:
        """
        Add a banner which shows the logo and summarise the number of FASTQ files
        processed

        """

        # Create the component
        self.expt_summary = ThroughputSummary(
            summary_name=self.summary_name,
            component_id="thoughput-summary",
            throughput_csv=throughput_csv,
        )

        # Define banner layout
        banner = html.Div(className="banner", children=self.expt_summary.get_layout())

        # Add to components and layout
        self.components.append(self.expt_summary)
        self.layout.append(banner)

    def _add_samples(self, samples_csv: str, samples_amplicons_csv: str) -> None:
        """
        Add a panel that shows progress of samples

        """
        self.samples = SamplesPie(
            summary_name=self.summary_name,
            component_id="samples-pie",
            samples_csv=samples_csv,
        )
        self.amplicons = AmpliconsBarplot(
            summary_name=self.summary_name,
            component_id="samples-barplot",
            samples_amplicons_csv=samples_amplicons_csv,
        )
        quality_row = html.Div(
            className="samples-row",
            children=[
                html.H3("Samples Statistics", style=dict(marginTop="0px")),
                html.Div(
                    className="samples-plots",
                    children=[self.samples.get_layout(), self.amplicons.get_layout()],
                ),
            ],
        )

        # Add components and layout
        self.components.append(self.samples)
        self.layout.append(quality_row)

    def _add_experiment_qc(self, coverage_csv: str) -> None:
        """
        Add a panel that shows quality control results

        """
        dropdown = dcc.Dropdown(
            id="quality-dropdown",
            options=[
                {"label": t(option), "value": option, "title": t(f"{option}_tooltip")}
                for option in QualityControl.STATISTICS
            ],
            value=QualityControl.STATISTICS[1],
            style=dict(width="300px"),
        )

        self.quality_control = QualityControl(
            self.summary_name,
            component_id="quality-heat",
            dropdown_id="quality-dropdown",
            coverage_csv=coverage_csv,
        )

        quality_row = html.Div(
            className="quality-row",
            children=[
                html.H3("Experiment QC Statistics", style=dict(marginTop="0px")),
                dropdown,
                html.Div(
                    className="quality-plots",
                    children=[self.quality_control.get_layout()],
                ),
            ],
        )

        # Add components and layout
        self.components.append(self.quality_control)
        self.layout.append(quality_row)

    def _add_prevalence_row(self, analysis_csv: str, master_csv: str) -> None:
        """
        Add a panel that shows prevalence calls

        """
        radio = dcc.RadioItems(
            id="prevalence-radio",
            options=list(PrevalenceBarplot.GENE_SETS.keys()),
            value=list(PrevalenceBarplot.GENE_SETS.keys())[0],
            inputClassName="prevalence-radio-input",
            labelClassName="prevalence-radio-label",
        )
        radio_by = dcc.RadioItems(
            id="prevalence-radio-by",
            options=["All", "region", "year", "region_year"],
            value="All",
            inputClassName="prevalence-radio-input",
            labelClassName="prevalence-radio-label",
        )

        self.prevalence_bars = PrevalenceBarplot(
            self.summary_name,
            component_id="prevalence-bars",
            radio_id="prevalence-radio",
            radio_id_by="prevalence-radio-by",
            analysis_csv=analysis_csv,
            master_csv=master_csv,
        )

        prevalence_row = html.Div(
            className="prevalence-row",
            children=[
                html.H3("Prevalence", style=dict(marginTop="0px")),
                html.Div(
                    className="prevalence-radio-row",
                    children=[radio, radio_by],
                ),
                html.Div(
                    className="prevalence-plots",
                    children=[self.prevalence_bars.get_layout()],
                ),
            ],
        )

        # Add components and layout
        self.components.append(self.prevalence_bars)
        self.layout.append(prevalence_row)

    def _add_prevalence_by_region_row(self, prevalence_region_csv: str) -> None:
        """
        Add a panel that shows prevalence calls by region

        """
        dropdown = dcc.Dropdown(
            id="gene-dropdown",
            options=PrevalenceBarplot.GENE_SETS["Resistance"],
            value=PrevalenceBarplot.GENE_SETS["Resistance"][0],
            style=dict(width="300px"),
        )

        self.prevalence_heatmap = PrevalenceHeatmap(
            summary_name=self.summary_name,
            prevalence_region_csv=prevalence_region_csv,
            component_id="prevalence-heatmap",
            gene_dropdown_id="gene-dropdown",
        )
        prevalence_row = html.Div(
            className="prevalence-region-row",
            children=[
                html.H3("Prevalence by region", style=dict(marginTop="0px")),
                dropdown,
                html.Div(
                    className="prevalence-region-plots",
                    children=[self.prevalence_heatmap.get_layout()],
                ),
            ],
        )

        # Add components and layout
        self.components.append(self.prevalence_heatmap)
        self.layout.append(prevalence_row)


class BasicSummaryDashboard(SummaryDashboardBuilder):
    """
    Build a dashboard with a focus on mapping statistics

    """

    CSS_STYLE = ["assets/summary-style.css"]

    def __init__(
        self,
        summary_name: str,
        throughput_csv: str,
        samples_csv: str,
        samples_amplicons_csv: str,
        coverage_csv: str,
        analysis_csv: str,
        master_csv: str,
        prevalence_region_csv: str,
    ):
        """
        Initialise all of the dashboard components

        """

        super().__init__(summary_name, self.CSS_STYLE)
        self.throughput_csv = throughput_csv
        self.samples_csv = samples_csv
        self.samples_amplicons_csv = samples_amplicons_csv
        self.coverage_csv = coverage_csv
        self.analysis_csv = analysis_csv
        self.master_csv = master_csv
        self.prevalence_region_csv = prevalence_region_csv

    def _gen_layout(self):
        """
        Generate the layout

        """
        self._add_throughput_banner(self.throughput_csv)
        self._add_samples(self.samples_csv, self.samples_amplicons_csv)
        self._add_experiment_qc(self.coverage_csv)
        self._add_prevalence_row(self.analysis_csv, self.master_csv)
        self._add_prevalence_by_region_row(self.prevalence_region_csv)
        # self._add_mapping_row(self.read_mapping_csv)
        # self._add_region_coverage_row(self.region_coverage_csv, self.regions)
        # self._add_depth_row(self.depth_profiles_csv, self.regions)
        # self._add_footer()
