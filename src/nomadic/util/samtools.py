import os
import uuid
import json
import subprocess


def samtools_view(input_bam, args, output_bam):
    """
    Run `samtools view` on an `input_bam`

    params
        input_bam : str
            Path to bam file.
        args : str
            Arguments to pass to `samtools view`
        output_bam : str
            Path to output bam file.

    returns
        None

    """

    cmd = "samtools view %s %s -o %s" % (input_bam, args, output_bam)
    subprocess.run(cmd, check=True, shell=True)

    return None


def samtools_index(input_bam):
    """
    Run index BAM

    params
        input_bam : str
            BAM file to be indexed.

    returns
        None

    """

    cmd = "samtools index %s" % input_bam
    subprocess.run(cmd, shell=True, check=True)

    return None


def samtools_merge(bam_files, output_bam):
    """
    Merge a collection of BAM files `bam_files`
    and write as an output BAM `output_bam`

    params
        bam_files : list, str, shape (n_bams, )
            A list of paths to bam files that should
            be merged.
        output_bam : str
            Path to output BAM file.

    returns
        None

    """

    cmd = "samtools merge -f %s %s" % (output_bam, " ".join(bam_files))
    subprocess.run(cmd, shell=True, check=True)

    return None


def samtools_flagstats(input_bam: str, output_json: str) -> None:
    """
    Run `samtools flagstats` on a target `bam`, and write the output
    to a JSON file

    """

    # Sanity checks, because `subprocess` can be opaque
    if not os.path.exists(input_bam):
        raise FileNotFoundError(f"Cannot find {input_bam}.")

    if not input_bam.endswith(".bam"):
        raise ValueError(f"Input bam file must end with `.bam`, currently: {input_bam}")

    # Create a temporary JSON
    temp_json = f"{input_bam[:-4]}.temp.{str(uuid.uuid4())[:8]}.json"

    # Run
    cmd = f"samtools flagstats -O json {input_bam} > {temp_json}"
    subprocess.run(cmd, shell=True, check=True)

    # Clean and write to output
    orig_dt = json.load(open(temp_json, "r"))["QC-passed reads"]
    clean_dt = {
        "n_reads": orig_dt["primary"],
        "n_mapped": orig_dt["mapped"],
        "n_uniq_mapped": orig_dt["primary mapped"],
        "n_chim_mapped": orig_dt["supplementary"],
        "n_mult_mapped": orig_dt["secondary"],
        "n_unmapped": orig_dt["primary"] - orig_dt["mapped"],
    }
    json.dump(clean_dt, open(output_json, "w"))

    # Remove temporary JSON
    os.remove(temp_json)
