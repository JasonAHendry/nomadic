# Small script to simulate a sequencing run;
# .fastqs are moved into appropriate directory at a given interval
# 2023/06/01, J. Hendry

import os
import shutil
import random
import time
import pandas as pd

AWAITED_FASTQ_DIR = "example_data/awaited_fastqs"
METADATA_PATH = "example_data/metadata/sample_info.csv"
FASTQ_PASS_DIR = "example_data/minknow/fastq_pass"
WAIT_INTERVAL_SECS = 5


def main():
    """
    Simulate the generation of FASTQ files as would occurr
    during a nanopore sequencing run

    """

    metadata = pd.read_csv(METADATA_PATH)
    barcodes = metadata["barcode"]

    barcode_fastqs = {}
    for barcode in barcodes:
        barcode_dir = f"{AWAITED_FASTQ_DIR}/{barcode}"
        fastqs = os.listdir(barcode_dir)
        barcode_fastqs[barcode] = [
            f"{barcode_dir}/{fastq}" for fastq in fastqs if fastq.endswith(".fastq.gz")
        ]

    while barcode_fastqs:
        b = random.choice(list(barcode_fastqs.keys()))
        f = random.choice(barcode_fastqs[b])

        print(f"{f} has been sequenced...")

        dst_dir = f"{FASTQ_PASS_DIR}/{os.path.basename(b)}"
        dst = f"{dst_dir}/{os.path.basename(f)}"
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)

        shutil.copyfile(src=f, dst=dst)

        barcode_fastqs[b].remove(f)
        if not barcode_fastqs[b]:
            barcode_fastqs.pop(b)

        time.sleep(WAIT_INTERVAL_SECS)

    print(f"No more FASTQ files in {AWAITED_FASTQ_DIR} to simulate sequencing.")
    clean = input("Clean up FASTQs for next time? [Yes/No]: ")
    if clean == "Yes":
        for barcode in barcodes:
            shutil.rmtree(f"{FASTQ_PASS_DIR}/{barcode}")
    else:
        print(f"Leaving FASTQs in {FASTQ_PASS_DIR}.")


if __name__ == "__main__":
    main()
