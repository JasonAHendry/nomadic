import os
import urllib.request
from nomadic.util.gff import load_gff


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

    def download_fasta(self):
        if self.ref.fasta_path and not self.exists_locally(self.ref.fasta_path):
            print("Downloading...")
            self.produce_dir(self.ref.fasta_path)
            urllib.request.urlretrieve(
                url=self.ref.fasta_url, filename=self.ref.fasta_path
            )
            print("Done.")
            print("")
        else:
            print("Already downloaded.")

    def download_gff(self, standardise: bool=False):
        if self.ref.gff_path and not self.exists_locally(self.ref.gff_path):
            print("Downloading...")
            self.produce_dir(self.ref.gff_path)
            urllib.request.urlretrieve(url=self.ref.gff_url, filename=self.ref.gff_path)
            print("Done.")
        else:
            print("Already downloaded.")

        if standardise:
            self._standardise_gff(self.ref.gff_path)

    def _standardise_gff(self, gff_path: str):
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
        gff_df = load_gff(gff_path)
        gff_df.query("feature in @KEEP_FIELDS", inplace=True)
        gff_df["feature"] = [
            to_gff3[f] if f in to_gff3 else f
            for f in gff_df["feature"]
        ]

        # Write to 'standardised' path
        gff_df.to_csv(self.ref.gff_standard_path, sep="\t", index=False, header=False)

