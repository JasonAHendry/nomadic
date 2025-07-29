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

    subprocess.run(["samtools", "view", input_bam, args, "-o", output_bam], check=True)

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

    subprocess.run(["samtools", "index", input_bam], check=True)

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

    subprocess.run(["samtools", "merge", "-f", output_bam, *bam_files], check=True)

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
    cmd = ["samtools", "flagstats", "-O", "json", input_bam]
    with open(temp_json, "w") as file:
        subprocess.run(cmd, check=True, stdout=file)

    # Clean and write to output
    orig_dt = json.load(open(temp_json, "r"))["QC-passed reads"]
    # NB: These are counting ALIGNMENTS not READS
    # A single read can have multiple alignments
    clean_dt = {
        "n_total": orig_dt["total"],
        "n_mapped": orig_dt["mapped"],
        "n_primary": orig_dt["primary mapped"],
        "n_secondary": orig_dt["secondary"],
        "n_supplementary": orig_dt["supplementary"],
        "n_unmapped": orig_dt["total"] - orig_dt["mapped"],
    }
    json.dump(clean_dt, open(output_json, "w"))

    # Remove temporary JSON
    os.remove(temp_json)


def samtools_depth(input_bam, output_path, region_str=None):
    """
    Run `samtools depth` on a given `input_bam`, focussing
    on regions defined by a `bed_file`

    params
        input_bam : str
            Path to bam file.
        output_path : str
            Path to write output file.
        region_bed : str
            BED file defining regions over which depth
            should be calculated. [optional]

    returns
        None

    """

    cmd = ["samtools", "depth"]
    if region_str is not None:
        cmd.extend(["-r", region_str])
    cmd.append("-aa")  # output all positions
    cmd.append("-J")
    cmd.extend(["-o", output_path])
    cmd.append(input_bam)
    subprocess.run(cmd, check=True)

    return None
