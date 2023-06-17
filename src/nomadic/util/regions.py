import pandas as pd
import seaborn as sns

from matplotlib.colors import rgb2hex

from nomadic.util.bed import load_bed_as_dataframe
from nomadic.util.exceptions import BEDFormatError


class RegionBEDParser:
    """
    Parse a BED file containing information about regions of interest,
    typically amplicons in targeted sequnencing

    """

    PALETTE = "tab20"

    def __init__(self, bed_path: str):
        """Load the BED file, assign colors to regions"""

        self.path = bed_path
        self.df = load_bed_as_dataframe(bed_path)

        self.n_regions = self.df.shape[0]
        self.names = self.df["name"].tolist()

        self._set_colors()

    def _check_names_unique(self):
        """Check all of the region names are unique"""

        if any(self.df["name"].duplicated()):
            raise BEDFormatError(
                "All regions must have unique names, but found duplicated."
            )

    def _set_colors(self):
        """
        Set the color scheme

        """

        self.colors = sns.color_palette(self.PALETTE, self.n_regions)
        self.colors_hex = [rgb2hex(c) for c in self.colors]
        self.col_map = dict(zip(self.names, self.colors))
        self.col_map_hex = dict(zip(self.names, self.colors_hex))
