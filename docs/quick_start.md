
To run *Nomadic* for a sequencing run, open a terminal window and type the following:

```
conda activate nomadic
cd <path/to/your/workspace>
nomadic realtime <expt_name>
```

- `<path/to/your/workspace>` should be replaced with the path to your *Nomadic* workspace.
- `<expt_name>` should be replaced with the name of your experiment.
    - You should have given your experiment the same name in *MinKNOW*.
    - You should have given your metadata file this name, and put it in your workspace metadata folder (`<path/to/your/workspace>/metadata/<expt_name>.csv|xlsx`).

The dashboard will open in a browser window on your computer. 

<!-- 1. Create a workspace if you have not already (`nomadic start pfalciparum`).
2. Pick an experiment name (e.g. `2025-06-21_first-sequencing`)
3. Start your sequencing run with *MinKNOW* using your experiment name.
4. Put your metadata file in the `metadata` folder in your *Nomadic* workspace. Name the file with your experiment name (e.g. `2025-06-21_first-sequencing.csv`).
5. Start *Nomadic*: `nomadic realtime 2025-06-21_first-sequencing` -->