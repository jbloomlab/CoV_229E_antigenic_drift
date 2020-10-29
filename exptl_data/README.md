# Analysis of 229E neutralization assay data

Each neutralization assay should have its data stored in an Excel file in a subdirectory named with the data in the format `2020-10-02`.
The subdirectories should also contain a `sample_map.csv` file that maps the Excel file data to the samples in a format that is readable by the Python script [excel_to_fractinfect.py](excel_to_fractinfect.py) (see [here](https://github.com/jbloomlab/exceltofractinfect) for more on this script, which was written by Kate Crawford).
The plate layouts referred to by the sample maps are in [./PlateLayouts/](PlateLayouts).

Then [Snakefile](Snakefile) can be used to analyze the neutralization data by running.
First make sure you have activated the `conda` environment with:

    conda activate CoV_229E_antigenic_drift

as described in the [top-level README](../README.md).
Then run [Snakefile](Snakefile) with the command:

    snakemake -j 1
