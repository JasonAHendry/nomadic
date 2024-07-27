import os
import urllib.request
from nomadic.util.gff import load_gff
from nomadic.util.fasta import find_lowcomplexity_intervals


class ReferenceDownloader:
    def __init__(self):
        self.ref = None

    def set_reference(self, reference):
        self.ref = reference

    @staticmethod
    def exists_locally(file_path):
        return os.path.isfile(file_path)

    @staticmethod
    def produce_dir(file_path):
        file_dir = os.path.dirname(file_path)
        if not os.path.isdir(file_dir):
            os.makedirs(file_dir)

    def download_fasta(self, create_mask: bool=False):
        if self.ref.fasta_path and not self.exists_locally(self.ref.fasta_path):
            print("Downloading FASTA...")
            print(f"  From: {self.ref.fasta_url}")
            print(f"  To: {self.ref.fasta_path}")
            self.produce_dir(self.ref.fasta_path)
            urllib.request.urlretrieve(
                url=self.ref.fasta_url, filename=self.ref.fasta_path
            )
            print("Done.")
            print("")
        else:
            print("Already downloaded FASTA.")

        if create_mask:
            if self.exists_locally(self.ref.fasta_mask_path):
                print("Already masked FASTA.")
            else:
                self._create_lowcomplexity_fasta_mask()

    def download_gff(self, standardise: bool=False):
        if self.ref.gff_path and not self.exists_locally(self.ref.gff_path):
            print("Downloading GFF...")
            print(f"  From: {self.ref.gff_url}")
            print(f"  To: {self.ref.gff_path}")
            self.produce_dir(self.ref.gff_path)
            urllib.request.urlretrieve(url=self.ref.gff_url, filename=self.ref.gff_path)
            print("Done.")
        else:
            print("Already downloaded GFF.")

        if standardise:
            if self.exists_locally(self.ref.gff_standard_path):
                print("Already standardised GFF.")
            else:
                self._standardise_gff()

    def _standardise_gff(self) -> None:
        """
        Try to standardise the GFF file into GFF3 format

        """
        
        # Settings
        KEEP_FIELDS = ["protein_coding_gene", "mRNA", "exon", "CDS"]
        to_gff3 = {
            "protein_coding_gene": "gene",
            "mRNA": "transcript"
        }

        # Standardise
        gff_df = load_gff(self.ref.gff_path)
        gff_df.query("feature in @KEEP_FIELDS", inplace=True)
        gff_df["feature"] = [
            to_gff3[f] if f in to_gff3 else f
            for f in gff_df["feature"]
        ]

        # Write to 'standardised' path
        gff_df.to_csv(self.ref.gff_standard_path, sep="\t", index=False, header=False)

    def _create_lowcomplexity_fasta_mask(self) -> None:
        """
        Create a BED file indicating regions that should be masked due to
        low complexity sequence
        
        """
        print("Creating a low-complexity mask for this reference genome (please be patient, this may take a few minutes)...")
        find_lowcomplexity_intervals(
            fasta_path=self.ref.fasta_path,
            bed_path=self.ref.fasta_mask_path
        )
        print("Done.\n")
