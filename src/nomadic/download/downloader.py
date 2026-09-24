import gzip
import os
import shutil
import time
import urllib.error
import urllib.request

from nomadic.download.references import Reference
from nomadic.util.fasta import find_lowcomplexity_intervals

# Retryable network/IO errors; a corrupt/truncated gzip stream also warrants a retry
RETRYABLE_ERRORS = (
    urllib.error.URLError,
    TimeoutError,
    ConnectionError,
    OSError,
    EOFError,
    gzip.BadGzipFile,
)


class ReferenceDownloader:
    MAX_RETRIES = 3
    RETRY_DELAY_SECS = 5

    def __init__(self):
        self.ref = None

    def set_reference(self, reference: Reference):
        self.ref = reference

    @staticmethod
    def exists_locally(file_path):
        return os.path.isfile(file_path)

    @staticmethod
    def produce_dir(file_path):
        file_dir = os.path.dirname(file_path)
        if not os.path.isdir(file_dir):
            os.makedirs(file_dir)

    def _download_with_retries(self, url: str, dest_path: str) -> None:
        """
        Download `url` to `dest_path`, retrying on transient failures and
        decompressing on the fly if `url` is gzipped. Downloads land in
        `.tmp` files named after `dest_path`, so `dest_path` itself is only
        ever created once the download (and any decompression) has fully
        succeeded.

        """
        self.produce_dir(dest_path)
        download_tmp_path = dest_path + ".download.tmp"
        decompressed_tmp_path = dest_path + ".tmp"
        needs_unzip = url.endswith(".gz") and not dest_path.endswith(".gz")

        last_error: BaseException = RuntimeError(
            f"Failed to download {url} after {self.MAX_RETRIES} attempts."
        )
        for attempt in range(1, self.MAX_RETRIES + 1):
            try:
                urllib.request.urlretrieve(
                    url=url, filename=download_tmp_path, reporthook=print_progress
                )
                print()  # Newline after progress bar which uses \r
                if needs_unzip:
                    print("Decompressing gzipped file...")
                    with (
                        gzip.open(download_tmp_path, "rb") as f_in,
                        open(decompressed_tmp_path, "wb") as f_out,
                    ):
                        shutil.copyfileobj(f_in, f_out)
                    os.remove(download_tmp_path)
                    os.replace(decompressed_tmp_path, dest_path)
                else:
                    os.replace(download_tmp_path, dest_path)
                return
            except RETRYABLE_ERRORS as e:
                last_error = e
                print(f"\nDownload attempt {attempt}/{self.MAX_RETRIES} failed: {e}")
                for tmp_path in (download_tmp_path, decompressed_tmp_path):
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
                if attempt < self.MAX_RETRIES:
                    time.sleep(self.RETRY_DELAY_SECS)
            except BaseException:
                for tmp_path in (download_tmp_path, decompressed_tmp_path):
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
                raise

        raise last_error

    def download_fasta(self, create_mask: bool = False):
        if self.ref is None:
            raise ValueError("Reference genome is not set.")
        if self.ref.fasta_path and not self.exists_locally(self.ref.fasta_path):
            print("Downloading FASTA...")
            print(f"  From: {self.ref.fasta_url}")
            print(f"  To: {self.ref.fasta_path}")
            self._download_with_retries(
                url=self.ref.fasta_url, dest_path=self.ref.fasta_path
            )
            print("Done.")
            print()
        else:
            print("Already downloaded FASTA.")

        if create_mask:
            if self.exists_locally(self.ref.fasta_mask_path):
                print("Already masked FASTA.")
            else:
                self._create_lowcomplexity_fasta_mask()

    def download_gff(self):
        if self.ref is None:
            raise ValueError("Reference genome is not set.")
        if self.ref.gff_path and not self.exists_locally(self.ref.gff_path):
            print("Downloading GFF...")
            print(f"  From: {self.ref.gff_url}")
            print(f"  To: {self.ref.gff_path}")
            self._download_with_retries(
                url=self.ref.gff_url, dest_path=self.ref.gff_path
            )
            print("Done.")
        else:
            print("Already downloaded GFF.")

    def _create_lowcomplexity_fasta_mask(self) -> None:
        """
        Create a BED file indicating regions that should be masked due to
        low complexity sequence

        """
        print(
            f"Creating a low-complexity mask at {self.ref.fasta_mask_path} for this reference genome (please be patient, this may take a few minutes)..."
        )
        tmp_path = self.ref.fasta_mask_path + ".tmp"
        find_lowcomplexity_intervals(fasta_path=self.ref.fasta_path, bed_path=tmp_path)
        os.rename(tmp_path, self.ref.fasta_mask_path)
        print("Done.\n")


def print_progress(block_number: int, block_size: int, total_size: int):
    downloaded = block_number * block_size
    if total_size > 0:
        percent = min(100, downloaded * 100 / total_size)
        print(f"Download progress: {downloaded} bytes ({percent:.2f}%)", end="\r")
    else:
        print(f"Download progress: {downloaded} bytes", end="\r")
