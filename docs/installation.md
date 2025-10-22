### Compatibility

_Nomadic_ runs on macOS or Linux. Windows is currently not supported.

### Requirements

To install _Nomadic_, you will need to install the package manager [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). _Mamba_ is faster and is recommended.

To install _Mamba_, run:

```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

or

```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

Then run the script with:

```
bash Miniforge3-$(uname)-$(uname -m).sh
```

### Steps

Open a terminal window and type:

```
conda create -c bioconda -n nomadic nomadic
```

Or, if you are using _Mamba_:

```
mambda create -c bioconda -n nomadic nomadic
```

_Nomadic_ should now be installed in it's own conda environment.

### Test installation

To test the installation, we will enter the conda environment:

```
conda activate nomadic
```

Then type:

```
nomadic --help
```

If `nomadic` has been installed successfully, you should see a set of available commands:

```
Usage: nomadic [OPTIONS] COMMAND [ARGS]...

  Mobile sequencing and analysis in real-time

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  start      Start a workspace.
  download   Download reference genomes.
  realtime   Run analysis in real-time.
  dashboard  Just run the dashboard.
```

### Update nomadic

To update _Nomadic_ to the newest version, activate your conda environment:

```
conda activate nomadic
```

Then type:

```
conda update nomadic
```

to update _Nomadic_.

You can check your installed version with

```
nomadic --version
```
