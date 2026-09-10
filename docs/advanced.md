## Realtime analysis

`nomadic realtime` analyses data while MinKNOW is still producing FASTQ files. It watches the FASTQ directory for each barcode in the experiment metadata, processes newly available files, and updates the experiment summaries. The command keeps watching until it is stopped with `Ctrl+C`.

The basic workflow is:

1. MinKNOW writes FASTQ files into the experiment's `fastq_pass` directory.
2. *Nomadic* finds new files for each barcode and maps them to the selected reference.
3. Mapping statistics, amplicon coverage, and depth profiles are updated.
4. If a caller is selected, SNP calls are also updated.
5. The dashboard and experiment-level summary files are refreshed after barcode updates.

The command runs a mapping and quality-control pipeline by default. Add `--caller bcftools` or `--caller delve` to include variant calling.

## Before you start

Create a workspace with `nomadic start`, then provide:

- A metadata CSV or XLSX file containing at least `barcode` and `sample_id` columns.
- A BED file describing the amplicons or regions of interest.
- A downloaded reference genome supported by *Nomadic*.
- A MinKNOW experiment whose name exactly matches the experiment name used by the command.

Inside a workspace, the usual invocation is:

```
nomadic realtime <experiment_name>
```

The metadata file is normally found at `<workspace>/metadata/<experiment_name>.csv`. The output is written to `<workspace>/results/<experiment_name>`.

## Custom panels and paths

Use `--region_bed` with either the name of a BED file in the workspace `beds` directory or a path to a custom BED file:

```
nomadic realtime experiment-01 \
    --region_bed path/to/my-panel.bed \
    --reference_name Pf3D7
```

You can provide every input explicitly when the files are outside the usual workspace layout:

```
nomadic realtime experiment-01 \
    --workspace /data/nomadic \
    --metadata_path /data/metadata/experiment-01.csv \
    --minknow_dir /var/lib/minknow/data \
    --output /data/results/experiment-01 \
    --region_bed /data/panels/my-panel.bed \
    --reference_name Pf3D7
```

When `--fastq_dir` is provided, it is used instead of `--minknow_dir`. Prefer `--minknow_dir` when possible because it allows *Nomadic* to retain the MinKNOW directory information used by related workflows.

## Workspace configuration

Each workspace can contain a `.config.yaml` file. *Nomadic* loads this file when the workspace is selected, either because the command is run inside the workspace or because `--workspace` points to it. Configuration values provide defaults; an option written on the command line always takes precedence.

The configuration file can define workspace-wide defaults and realtime-specific defaults:

```yaml
defaults:
  region_bed: nomadsMVP
  reference_name: Pf3D7
  caller: delve
  minknow_dir: /var/lib/minknow/data

realtime:
  defaults:
    threads: 8
    dashboard: true
```

Configuration keys are created from the argument name of the command-line option. Use the long option name without its leading `--`, replacing hyphens (`-`) with underscores (`_`). Short options such as `-b` are not used as keys. For example:

| Command-line option | Configuration key |
| --- | --- |
| `--region_bed` | `region_bed` |
| `--reference-name` | `reference_name` |
| `--minknow_dir` | `minknow_dir` |
| `--no-dashboard` | `dashboard` |

The key must be placed in the defaults section that should provide it: `defaults` for a workspace-wide default or `realtime.defaults` for a realtime-only default.

The top-level `defaults` section is shared by commands. Values under `realtime.defaults` apply only to `nomadic realtime`. When both sections define the same option, the realtime-specific value is used. Explicit command-line options override both sections:

```
nomadic realtime experiment-01 --caller bcftools
```

In this example, `bcftools` is used even though the configuration file specifies `delve`.

When a workspace is created with `nomadic start`, its `.config.yaml` is initialized with the organism's default `region_bed`, `reference_name`, and `caller`. Edit that file to change those defaults or add values such as `minknow_dir`. The file is ordinary YAML and can be edited directly.

!!! note
    A workspace is needed when an option depends on workspace defaults. If the current directory is not a workspace, provide `--workspace` or specify the required paths explicitly.

## Command options

### Input and output

| Option | Default | Description |
| --- | --- | --- |
| `EXPERIMENT_NAME` | Required | Name of the MinKNOW experiment. It must match the experiment name used when sequencing and the metadata filename. |
| `-w`, `--workspace` | Current directory | Workspace containing the `beds`, `metadata`, and `results` directories. |
| `-m`, `--metadata_path` | `<workspace>/metadata/<experiment_name>.csv` | Metadata CSV or XLSX containing barcode and sample information. |
| `-o`, `--output` | `<workspace>/results/<experiment_name>` | Directory where this experiment's results are stored. |
| `-k`, `--minknow_dir` | MinKNOW's default data directory | MinKNOW base directory or experiment directory. |
| `-f`, `--fastq_dir` | Resolved from `--minknow_dir` | FASTQ directory or glob. If provided, it takes precedence over `--minknow_dir`. |

### Reference and regions

| Option | Default | Description |
| --- | --- | --- |
| `-b`, `--region_bed` | Required unless configured | BED file path or panel name, such as `nomads8` or `nomadsMVP`. |
| `-r`, `--reference_name` | Required unless configured | Reference genome used for mapping and downstream analysis. Supported values include `Pf3D7`, `PfDd2`, `Pv`, `Poc`, `Pm`, `AgPEST`, `AaDONGOLA2021`, `AcolN3`, `AfunGA1`, `AsUCISS2018`, and `Hs`. |

The selected reference must already be available to *Nomadic*. Reference availability is checked before processing begins.

### Variant calling

| Option | Default | Description |
| --- | --- | --- |
| `-c`, `--caller` | No variant calling | Select `bcftools` or `delve` for SNP calling. |

Variant calling adds VCF and variant summary outputs to the mapping and coverage results. The detailed output columns are described in [Output files](output_files.md).

### Run control

| Option | Default | Description |
| --- | --- | --- |
| `--resume` | `false` | Explicitly resume an existing experiment output directory. Normally *Nomadic* prompts when the output directory already exists. |
| `--overwrite` | `false` | Delete the existing output directory and start again. All results in that directory are removed. |
| `-v`, `--verbose` | `false` | Increase logging verbosity for debugging. |

If a run is interrupted, start the same command again and choose to resume, or pass `--resume`. *Nomadic* records completed FASTQ increments in a `.work.log` file for each barcode and reprocesses any increment that was started but not completed. It also checks that the inputs match the original run before continuing.

### Performance

| Option | Default | Description |
| --- | --- | --- |
| `-t`, `--threads` | `5` | Number of analysis threads. Increasing this can improve throughput, but also increases CPU and memory use. |

## Output and recovery

An experiment output directory contains experiment-level summary files, a `metadata` directory, and one directory for each barcode. The metadata directory includes a copy of the parsed metadata, the selected BED file, and `settings.json`, which records the inputs used for the run.

The summary files include:

- `summary.fastqs_processed.csv`
- `summary.read_mapping.csv`
- `summary.region_coverage.csv`
- `summary.depth_profiles.csv`
- `summary.aa_changes.csv` and `summary.nt_changes.csv` when variant calling is enabled

Results are updated incrementally as new FASTQ files arrive. The watcher records the start and completion of each increment in `<barcode>/.work.log`, allowing an interrupted increment to be detected and rerun after restart. See [Output files](output_files.md) for the columns and interpretation of the generated tables.

## Troubleshooting

### The workspace cannot be found

Run the command from inside a workspace, pass `--workspace <path>`, or provide explicit values for the metadata, BED, reference, and output options.

### No FASTQ files are processed

Check that the experiment name matches MinKNOW exactly and that the expected barcode directories are inside the resolved FASTQ directory. Use `--minknow_dir` to point to the MinKNOW base directory or `--fastq_dir` to provide the FASTQ location directly.

### The output directory already exists

Use `--resume` (or answer `y` when prompted) to force continuation of an existing run. Use `--overwrite` (or respond `r` when prompted) only when the previous results should be deleted and the experiment should be started again.

### The configuration is rejected

Check that `.config.yaml` is valid YAML, that `defaults` and `realtime.defaults` are mappings, and that option values use the expected types. Run with an explicit command-line option to temporarily override a configured value while diagnosing the file.