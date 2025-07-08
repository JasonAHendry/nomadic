from .downloader import ReferenceDownloader
from .references import REFERENCE_COLLECTION
from .corrections import DHPS, update_reference_genome


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
            downloader.download_fasta(create_mask=True)
            downloader.download_gff(standardise=True)

            if r.name == "Pf3D7":  # convert DHPS to WT
                update_reference_genome(r.fasta_path, [DHPS])

    elif reference_name is not None:
        print(f"Reference: {reference_name}")
        r = REFERENCE_COLLECTION[reference_name]
        downloader.set_reference(r)
        downloader.download_fasta(create_mask=True)
        downloader.download_gff(standardise=True)

        if r.name == "Pf3D7":  # convert DHPS to WT
            update_reference_genome(r.fasta_path, [DHPS])
    else:
        print("Must specify options. Type 'nomadic download --help' for details.")
