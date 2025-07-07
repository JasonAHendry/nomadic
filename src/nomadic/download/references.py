import os
from nomadic.util.exceptions import ReferenceGenomeMissingError
from abc import ABC, abstractmethod


# ================================================================
# Base classes for reference sequences
#
# ================================================================


class Reference(ABC):
    """
    Basic interface for reference sequences used for
    mapping, variant calling, &c.

    """

    def __init__(self):
        self.name = None
        self.fasta_url = None
        self.gff_url = None
        self.fasta_path = None
        self.gff_path = None

    @property
    def gff_standard_path(self):
        if self.gff_path is None:
            return None
        return self.gff_path.replace(".gff", ".standard.gff")

    @property
    def fasta_mask_path(self):
        if self.fasta_path is None:
            return None
        return self.fasta_path.replace(".fasta", ".lowcomplexity_mask.bed")

    @abstractmethod
    def set_fasta(self) -> None:
        pass

    @abstractmethod
    def set_gff(self) -> None:
        pass

    @staticmethod
    def exists_locally(file_path):
        return os.path.isfile(file_path)

    def confirm_downloaded(self) -> None:
        """
        Confirm that a reference genome has been downloaded,
        raising a ReferenceGenomeMissingError if not.

        Useful to invoke before running analysis steps that depend
        on the reference being present (e.g. mapping).

        """

        if not self.exists_locally(self.fasta_path):
            raise ReferenceGenomeMissingError(
                f"For the reference genome '{self.name}'"
                + f" the FASTA file is missing. Please run `nomadic download -r {self.name}`."
            )

        if not self.exists_locally(self.gff_path):
            raise ReferenceGenomeMissingError(
                f"For the reference genome '{self.name}'"
                + f" the GFF file is missing. Please run `nomadic download -r {self.name}`."
            )

        if not self.exists_locally(self.gff_standard_path):
            raise ReferenceGenomeMissingError(
                f"For the reference genome '{self.name}'"
                + f" the standardised GFF file is missing. Please run `nomadic download -r {self.name}`."
            )

        if not self.exists_locally(self.fasta_mask_path):
            raise ReferenceGenomeMissingError(
                f"For the reference genome '{self.name}'"
                + f" the masked FASTA file is missing. Please run `nomadic download -r {self.name}`."
            )


class PlasmoDB(Reference):
    """
    Encapsulate reference sequence downloads from PlasmoDB

    """

    source = "plasmodb"
    source_url = "https://plasmodb.org/common/downloads"
    release = 67

    def __init__(self, species, strain):
        self.species = species
        self.strain = strain
        self.data_url = f"{self.source_url}/release-{self.release}/{species}{strain}"
        self.set_fasta()
        self.set_gff()

    def set_fasta(self):
        """Set .fasta file download URL and local path"""
        fasta_fn = f"PlasmoDB-{self.release}_{self.species}{self.strain}_Genome.fasta"
        self.fasta_url = f"{self.data_url}/fasta/data/{fasta_fn}"
        self.fasta_path = f"resources/{self.source}/{self.release}/{fasta_fn}"

    def set_gff(self):
        """Set .gff file download URL and local path"""

        gff_fn = f"PlasmoDB-{self.release}_{self.species}{self.strain}.gff"
        self.gff_url = f"{self.data_url}/gff/data/{gff_fn}"
        self.gff_path = f"resources/{self.source}/{self.release}/{gff_fn}"


class VectorBase(Reference):
    """
    Encapsulate reference sequence downloads from vectorbase

    """

    source = "vectorbase"
    source_url = "https://vectorbase.org/common/downloads"
    release = 67

    def __init__(self, species, strain):
        self.species = species
        self.strain = strain
        self.data_url = f"{self.source_url}/release-{self.release}/{species}{strain}"
        self.set_fasta()
        self.set_gff()

    def set_fasta(self):
        """Set .fasta file download URL and local path"""
        fasta_fn = f"VectorBase-{self.release}_{self.species}{self.strain}_Genome.fasta"
        self.fasta_url = f"{self.data_url}/fasta/data/{fasta_fn}"
        self.fasta_path = f"resources/{self.source}/{self.release}/{fasta_fn}"

    def set_gff(self):
        """Set .gff file download URL and local path"""

        gff_fn = f"VectorBase-{self.release}_{self.species}{self.strain}.gff"
        self.gff_url = f"{self.data_url}/gff/data/{gff_fn}"
        self.gff_path = f"resources/{self.source}/{self.release}/{gff_fn}"


# ===============================================================
# Classes for specific reference sequences
#
# ================================================================


class PlasmodiumFalciparum3D7(PlasmoDB):
    def __init__(self):
        self.name = "Pf3D7"
        super().__init__(species="Pfalciparum", strain="3D7")


class PlasmodiumFalciparumDd2(PlasmoDB):
    def __init__(self):
        self.name = "PfDd2"
        super().__init__(species="Pfalciparum", strain="Dd2")


class PlasmodiumVivax(PlasmoDB):
    def __init__(self):
        self.name = "Pv"
        super().__init__(species="Pvivax", strain="P01")


class PlasmodiumOvale(PlasmoDB):
    def __init__(self):
        self.name = "Poc"
        super().__init__(species="Povalecurtisi", strain="GH01")


class PlasmodiumMalariae(PlasmoDB):
    def __init__(self):
        self.name = "Pm"
        super().__init__(species="Pmalariae", strain="UG01")


class AnophelesGambiaePEST(VectorBase):
    def __init__(self):
        self.name = "AgPEST"
        super().__init__(species="Agambiae", strain="PEST")


class AnophelesArabiensis(VectorBase):
    def __init__(self):
        self.name = "AaDONGOLA2021"
        super().__init__(species="Aarabiensis", strain="DONGOLA2021")


class AnophelesColuzzi(VectorBase):
    def __init__(self):
        self.name = "AcolN3"
        super().__init__(species="Acoluzzii", strain="AcolN3")


class AnophelesFunestus(VectorBase):
    def __init__(self):
        self.name = "AfunGA1"
        super().__init__(species="Afunestus", strain="AfunGA1")


class AnophelesStephensi(VectorBase):
    def __init__(self):
        self.name = "AsUCISS2018"
        super().__init__(species="Astephensi", strain="UCISS2018")


class HomoSapiens(Reference):
    """
    Download the Homo Sapiens reference genome from
    Ensembl

    """

    source = "ensembl"
    source_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
    source_url += "000/001/405/GCA_000001405.15_GRCh38/"

    def __init__(self):
        self.name = "Hs"
        self.set_fasta()
        self.set_gff()

    def set_fasta(self):
        fasta_fn = "GCA_000001405.29_GRCh38.p14_genomic.fna.gz"
        self.fasta_url = f"{self.source_url}/{fasta_fn}"
        self.fasta_path = f"resources/{self.source}/{fasta_fn}"

    def set_gff(self):
        gff_fn = "GCA_000001405.29_GRCh38.p14_genomic.gff.gz"
        self.gff_url = f"{self.source_url}/{gff_fn}"
        self.gff_path = f"resources/{self.source}/{gff_fn}"


# ================================================================
# Collection
# Note they are already initialised
# ================================================================


REFERENCE_COLLECTION = {
    r.name: r
    for r in [
        PlasmodiumFalciparum3D7(),
        PlasmodiumFalciparumDd2(),
        PlasmodiumVivax(),
        PlasmodiumOvale(),
        PlasmodiumMalariae(),
        AnophelesGambiaePEST(),
        AnophelesArabiensis(),
        AnophelesColuzzi(),
        AnophelesFunestus(),
        AnophelesStephensi(),
        HomoSapiens(),
    ]
}
