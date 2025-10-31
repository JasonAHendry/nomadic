import glob
from importlib.resources import as_file, files
import logging
import threading
from abc import ABC, abstractmethod
from dash import Dash, html, dcc

from i18n import t
import i18n
import pandas as pd

# from importlib.resources import files, as_file
from nomadic.summarize.dashboard.components import (
    AmpliconsBarplot,
    GeneDeletionsBarplot,
    PrevalenceHeatmap,
    SamplesPie,
    ThroughputSummary,
    QualityControl,
    PrevalenceBarplot,
    MapComponent,
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

        setup_translations()
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
                html.H3("Sample Statistics", style=dict(marginTop="0px")),
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
                html.Div(
                    className="quality-dropdowns",
                    children=[
                        html.Div(
                            children=[
                                html.Label("Select statistic:"),
                                dropdown,
                            ]
                        ),
                    ],
                ),
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
        dropdown_genset = dcc.Dropdown(
            id="prevalence-dropdown-gene-set",
            options=list(PrevalenceBarplot.GENE_SETS.keys()),
            value=list(PrevalenceBarplot.GENE_SETS.keys())[0],
            style=dict(width="300px"),
            clearable=False,
        )

        cols = cols_to_group_by(master_csv, analysis_csv, 10)

        dropdown_by = dcc.Dropdown(
            id="prevalence-dropdown-by",
            options=["All", *cols],
            value="All",
            style=dict(width="300px"),
            clearable=False,
        )

        self.prevalence_bars = PrevalenceBarplot(
            self.summary_name,
            component_id="prevalence-bars",
            radio_id="prevalence-dropdown-gene-set",
            radio_id_by="prevalence-dropdown-by",
            analysis_csv=analysis_csv,
            master_csv=master_csv,
        )

        prevalence_row = html.Div(
            className="prevalence-row",
            children=[
                html.H3("Prevalence", style=dict(marginTop="0px")),
                html.Div(
                    className="prevalence-dropdowns",
                    children=[
                        html.Div(
                            children=[
                                html.Label("Select gene set:"),
                                dropdown_genset,
                            ]
                        ),
                        html.Div(
                            children=[
                                html.Label("Group by:"),
                                dropdown_by,
                            ]
                        ),
                    ],
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

    def _add_prevalence_by_col_row(self, analysis_csv: str, master_csv: str) -> None:
        """
        Add a panel that shows prevalence calls by cols

        """
        gen_dropdown = dcc.Dropdown(
            id="gene-dropdown",
            options=PrevalenceBarplot.GENE_SETS["Resistance"],
            value=PrevalenceBarplot.GENE_SETS["Resistance"][0],
            style=dict(width="300px"),
            clearable=False,
        )

        cols = cols_to_group_by(master_csv, analysis_csv, 50)

        col_dropdown = dcc.Dropdown(
            id="col-dropdown",
            options=cols,
            value=cols[0],
            style=dict(width="300px"),
            clearable=False,
        )

        self.prevalence_heatmap = PrevalenceHeatmap(
            summary_name=self.summary_name,
            analysis_csv=analysis_csv,
            master_csv=master_csv,
            component_id="prevalence-heatmap",
            gene_dropdown_id="gene-dropdown",
            col_dropdown_id="col-dropdown",
        )
        prevalence_row = html.Div(
            className="prevalence-by-row",
            children=[
                html.H3("Prevalence by category", style=dict(marginTop="0px")),
                html.Div(
                    className="prevalence-dropdowns",
                    children=[
                        html.Div(
                            children=[
                                html.Label("Select gene:"),
                                gen_dropdown,
                            ]
                        ),
                        html.Div(
                            children=[
                                html.Label("Group by:"),
                                col_dropdown,
                            ]
                        ),
                    ],
                ),
                html.Div(
                    className="prevalence-region-plots",
                    children=[self.prevalence_heatmap.get_layout()],
                ),
            ],
        )

        # Add components and layout
        self.components.append(self.prevalence_heatmap)
        self.layout.append(prevalence_row)

    def _add_gene_deletion_row(self, gene_deletions_csv: str, master_csv: str) -> None:
        """
        Add a panel that shows prevalence calls

        """

        cols = cols_to_group_by(master_csv, gene_deletions_csv, 50)

        dropdown_by = dcc.Dropdown(
            id="gene-deletions-dropdown-by",
            options=["All", *cols],
            value="All",
            style=dict(width="300px"),
            clearable=False,
        )

        self.prevalence_bars = GeneDeletionsBarplot(
            self.summary_name,
            component_id="gene-deletions-bars",
            radio_id_by="gene-deletions-dropdown-by",
            gene_deletions_csv=gene_deletions_csv,
            master_csv=master_csv,
        )

        prevalence_row = html.Div(
            className="gene-deltions-row",
            children=[
                html.H3("Prevalence Gene Deletions", style=dict(marginTop="0px")),
                html.Div(
                    className="prevalence-dropdowns",
                    children=[
                        html.Div(
                            children=[
                                html.Label("Group by:"),
                                dropdown_by,
                            ]
                        ),
                    ],
                ),
                html.Div(
                    className="gene-deletions-plots",
                    children=[self.prevalence_bars.get_layout()],
                ),
            ],
        )

        # Add components and layout
        self.components.append(self.prevalence_bars)
        self.layout.append(prevalence_row)

    def _add_map_row(
        self, analysis_csv: str, master_csv: str, geojsons: list[str]
    ) -> None:
        """
        Add a panel that shows a choropleth map of drug resistance marker prevalence
        """
        # Get unique mutations from the analysis CSV for resistance genes only
        analysis_df = pd.read_csv(analysis_csv)
        resistance_genes = PrevalenceBarplot.GENE_SETS["Resistance"]
        resistance_df = analysis_df[analysis_df["gene"].isin(resistance_genes)]
        resistance_df["gene_mutation"] = (
            resistance_df["gene"] + "-" + resistance_df["aa_change"]
        )
        gene_mutations = sorted(resistance_df["gene_mutation"].unique())

        # Create the dropdowns
        mutation_dropdown = dcc.Dropdown(
            id="map-mutation-dropdown",
            options=gene_mutations,
            value=gene_mutations[0] if gene_mutations else None,
            style=dict(width="300px"),
            clearable=False,
        )

        regions = {
            path.split("/")[-1].split(".")[0].split("-")[1]: path for path in geojsons
        }

        region_dropdown = dcc.Dropdown(
            id="map-region-dropdown",
            options=list(regions.keys()),
            value="district",
            style=dict(width="300px"),
            clearable=False,
        )

        # Create the map component with the prepared dropdowns
        self.prevalence_map = MapComponent(
            summary_name=self.summary_name,
            analysis_csv=analysis_csv,
            master_csv=master_csv,
            component_id="prevalence-map",
            mutation_dropdown_id="map-mutation-dropdown",
            region_dropdown_id="map-region-dropdown",
            geojsons=regions,
        )

        map_row = html.Div(
            className="map-row",
            children=[
                html.H3("Geographic Distribution", style=dict(marginTop="0px")),
                html.Div(
                    className="map-dropdowns",
                    children=[
                        html.Div(
                            children=[
                                html.Label("Select mutation:"),
                                mutation_dropdown,
                            ]
                        ),
                        html.Div(
                            children=[
                                html.Label("Group by:"),
                                region_dropdown,
                            ]
                        ),
                    ],
                ),
                html.Div(
                    className="map-plot",
                    children=[self.prevalence_map.get_layout()],
                ),
            ],
        )

        # Add components and layout
        self.components.append(self.prevalence_map)
        self.layout.append(map_row)


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
        gene_deletions_csv: str,
        master_csv: str,
        geojson_glob: str,
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
        self.gene_deletions_csv = gene_deletions_csv
        self.geojson_glob = geojson_glob

    def _gen_layout(self):
        """
        Generate the layout

        """
        self._add_throughput_banner(self.throughput_csv)
        self._add_samples(self.samples_csv, self.samples_amplicons_csv)
        self._add_experiment_qc(self.coverage_csv)
        self._add_prevalence_row(self.analysis_csv, self.master_csv)
        self._add_prevalence_by_col_row(self.analysis_csv, self.master_csv)
        self._add_gene_deletion_row(self.gene_deletions_csv, self.master_csv)

        if glob.glob(self.geojson_glob):
            self._add_map_row(
                self.analysis_csv, self.master_csv, glob.glob(self.geojson_glob)
            )


def setup_translations():
    """
    Set up translations for the dashboard

    This function loads the translation files from the package resources
    and appends them to the i18n load path.

    """
    with as_file(files("nomadic.summarize.dashboard").joinpath("translations")) as path:
        i18n.load_path.append(str(path))
        i18n.set("filename_format", "{locale}.{format}")
        i18n.load_everything()
    i18n.set("locale", "en")


def cols_to_group_by(master_csv: str, analysis_csv, max_cat: int) -> list[str]:
    """
    Get columns that can be used to group prevalence by

    """

    master_df = pd.read_csv(master_csv)
    analysis_df = pd.read_csv(analysis_csv)
    df = pd.merge(
        master_df,
        analysis_df[["sample_id"]],
        on="sample_id",
        how="inner",
    )
    cols = df.columns.tolist()
    cols.remove("sample_id")

    for col in cols[:]:
        if pd.api.types.is_numeric_dtype(df[col]):
            cols.remove(col)
            continue
        n_unique = df[col].nunique()
        if n_unique > max_cat:
            cols.remove(col)

    return cols
