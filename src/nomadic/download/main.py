from nomadic.util.logging_config import LoggingFascade
from .downloader import ReferenceDownloader
from .references import REFERENCE_COLLECTION


def main(reference_name: str, all: bool = False) -> None:
    """
    Main script for downloading various reference genomes defined within
    the `REFERENCE_COLLECTION`

    """

    downloader = ReferenceDownloader()

    if all:
        print(
            f"Will download all {len(REFERENCE_COLLECTION)} available reference genomes."
        )
        for rname, r in REFERENCE_COLLECTION.items():
            print(f"Reference: {rname}")
            downloader.set_reference(r)
            downloader.download_fasta()
            downloader.download_gff(standardise=True) # TODO: Might not want always
    elif reference_name is not None:
        print(f"Reference: {reference_name}")
        downloader.set_reference(REFERENCE_COLLECTION[reference_name])
        downloader.download_fasta()
        downloader.download_gff(standardise=True)
    else:
        print("Must specify options. Type 'nomadic download --help' for details.")
