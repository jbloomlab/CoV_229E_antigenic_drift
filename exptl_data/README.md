# Analysis of titer neutralization assay data
This subdirectory contains data and analysis relevant to the experiments.

The maps of the spike-expressing plasmids are in [./plasmid_maps/](plasmid_maps).

There is also code that analyzes and fits the neutralization data.
Each neutralization assay should have its data stored in an Excel file in a subdirectory named with the data in the format `2020-10-02`.
The subdirectories should also contain a `sample_map.csv` file that maps the Excel file data to the samples in a format that is readable by the Python script [excel_to_fractinfect.py](excel_to_fractinfect.py) (see [here](https://github.com/jbloomlab/exceltofractinfect) for more on this script, which was written by Kate Crawford).
The plate layouts referred to by the sample maps are in [./PlateLayouts/](PlateLayouts).

The file [serum_info.csv](serum_info.csv) maps the serum names to relevant information about the individual.

Then [Snakefile](Snakefile) can be used to analyze the neutralization data.
It does this by first processing the Excel neutralization data, and then running the Jupyter notebook [analyze_neut_data.ipynb](analyze_neut_data.ipynb).
First make sure you have activated the `conda` environment with:

    conda activate CoV_229E_antigenic_drift

as described in the [top-level README](../README.md).
Then run [Snakefile](Snakefile) with the command:

    snakemake -j 1

The results of the analysis are in the Jupyter notebook [analyze_neut_data.ipynb](analyze_neut_data.ipynb) and its Markdown rendering at [results/analyze_neut_data.md](results/analyze_neut_data.md).

We also plot the viral titers in [viral_titers.csv](viral_titers.csv) using the notebook [plot_titers.ipynb](plot_titers.ipynb), which is rendered in Markdown at [results/plot_titers.md](results/plot_titers.md)
