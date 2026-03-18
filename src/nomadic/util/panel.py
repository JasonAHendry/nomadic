from dataclasses import dataclass


@dataclass
class PanelSettings:
    """Settings for a panel in the summary."""

    name: str
    # List of amplicons to exclude from analysis set
    excluded_amplicons: list[str]
    # List of mutations to exclude from analysis set
    filtered_mutations: list[str]
    # Amplicon sets
    amplicon_sets: dict[str, list[str]]
    # List of genes for which to compute deletion prevalence
    deletion_genes: list[str]


MVP_PANEL_SETTINGS = PanelSettings(
    name="nomadsMVP",
    excluded_amplicons=[
        "hrp2-p14-306",
        "hrp3-p14-276",
    ],
    filtered_mutations=["crt-N75K"],
    amplicon_sets={
        "Resistance": [
            "crt-p14-125",
            "dhps-p317-707",
            "dhfr-p1-410",
            "kelch13-p383-727",
            "mdr1-p968-1278",
            "mdr1-p46-245",
        ],
        "Diversity": ["ama1-p74-384", "csp-p19-398"],
    },
    # Don't show deletion genes at the moment, because the code is not robust
    # deletion_genes=["hrp2", "hrp3"],
    deletion_genes=[],
)


UNKNOWN_PANEL_SETTINGS = PanelSettings(
    name="Unknown",
    excluded_amplicons=[],
    filtered_mutations=[],
    amplicon_sets={},
    deletion_genes=[],
)


def get_panel_settings(panel_name: str) -> PanelSettings:
    """Get panel settings by panel name."""
    if panel_name == "nomadsMVP":
        return MVP_PANEL_SETTINGS
    else:
        return UNKNOWN_PANEL_SETTINGS
