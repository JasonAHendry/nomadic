import os
import urllib.request


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

    def download_gff(self):
        if self.ref.gff_path and not self.exists_locally(self.ref.gff_path):
            print("Downloading...")
            self.produce_dir(self.ref.gff_path)
            urllib.request.urlretrieve(url=self.ref.gff_url, filename=self.ref.gff_path)
            print("Done.")
        else:
            print("Already downloaded.")
