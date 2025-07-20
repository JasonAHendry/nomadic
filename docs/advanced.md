We will populate this page with more details soon. But for now, we would direct advanced users to the [*Nomadic* github](https://github.com/JasonAHendry/nomadic) and the help text of individual commands for advanced information about *Nomadic*.

For a quick refernece, here is the help for `nomadic realtime`

```
Usage: nomadic realtime [OPTIONS] EXPERIMENT_NAME

  Analyse data being produced by MinKNOW while sequencing is ongoing

Options:
  -w, --workspace DIRECTORY       Path of the workspace where all input/output
                                  files (beds, metadata, results) are stored.
                                  The workspace directory simplifies the use
                                  of nomadic in that many arguments don't need
                                  to be listed as they are predefined in the
                                  workspace config or can be loaded from the
                                  workspace  [default: (current directory)]
   TEXT
  -o, --output PATH               Path to the output directory where results
                                  of this experiment will be stored. Usually
                                  the default of storing it in the workspace
                                  should be enough.  [default:
                                  (<workspace>/results/<experiment_name>)]
  -f, --fastq_dir DIRECTORY       Path to `fastq_pass` directory produced by
                                  MinKNOW or Guppy.  [default: (/var/lib/minkn
                                  ow/data/<experiment_name>/fastq_pass)]
  -m, --metadata_csv FILE         Path to metadata CSV file containing barcode
                                  and sample information.  [default: (<workspa
                                  ce>/metadata/<experiment_name>.csv)]
  -b, --region_bed PATH           Path to BED file specifying genomic regions
                                  of interest or name of panel, e.g. 'nomads8'
                                  or 'nomadsMVP'.  [required]
  -r, --reference_name [Pf3D7|PfDd2|Pv|Poc|Pm|AgPEST|AaDONGOLA2021|AcolN3|AfunGA1|AsUCISS2018|Hs]
                                  Choose a reference genome to be used in
                                  real-time analysis.
  -c, --call                      Perform preliminary variant calling of
                                  biallelic SNPs in real-time.
  --resume                        Resume a previous experiment run if the
                                  output directory already exists. Only use if
                                  you want to force resuming an already
                                  started experiment. Not needed in
                                  interactive mode as this will be prompted
  -v, --verbose                   Increase logging verbosity. Helpful for
                                  debugging.
  --help                          Show this message and exit.
```

You can see that almost all inputs to the software can be specified explicitly (instead of using the default behaviour). For example, to use your own amplicon panel with *P. falciparum* malaria, you just need to specify the path to it using the `-b` flag.