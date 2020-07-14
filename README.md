# Antigenic drift of the human seasonal coronavirus 229E

Study by Laurel Kelnhofer-Millevolte, Rachel Eguia, Katharine Crawford, and Jesse Bloom

## Steps in analysis

Most of the analysis is run automatically via [snakemake](https://snakemake.readthedocs.io/), but there is an initial manual step to get Genbank accessions and metadata.

### Manual step: get Genbank accessions and extra metadata
For the analysis, we need Genbank accessions and relevant meta data.
We download all relevant accessions from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/).
In addition, we manually specify metadata that is described in publications for some accessions that don't have the relevant metadata in the download from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/).

To download the accessions, go to the [NCBI Virus website](https://www.ncbi.nlm.nih.gov/labs/virus/), click *Search by virus*, and then enter *Human coronavirus 229E* (taxid:1137) in the *Search by virus name or taxonomy* box and hit return.
Then click *Download*, select *CSV format* under *Current table view results*, hit *Next*, hit *Next* again (it should show *Download all records*), and enter choose *Select All* for the columns to include.
Then save the file as [data/NCBI_Virus_229E_accessions.csv](data/NCBI_Virus_229E_accessions.csv).
The current file was downloaded on July 13, 2020.

Some accessions are missing information we need, such as about host species and collection date.
The file [data/extra_229E_accessions_metadata.yaml](data/extra_229E_accessions_metadata.yaml) has this missing data manually entered as parsed from publications.
The format of this file should be self-explanatory, and describes the sources of the manually parsed data.

### Steps run automatically by Snakemake pipeline
The remaining steps are run automatically using [snakemake](https://snakemake.readthedocs.io/) to run [Snakefile](Snakefile), which reads its configuration from [config.yaml](config.yaml).
The results are placed in [./results/](results).

To run these steps, first build the conda environment, which installs the necessary programs.
First [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
Then build the environment.
If you are using a Mac OS X, you can probably install the full version-pinned environment in [environment_pinned.yml](environment_pinned.yml) with:

    conda env create -f environment_pinned.yml

If this fails (which is especially likely if you are using a different operating system), then install the version-unpinned environment in [environment_unpinned.yml](environment_unpinned.yml), which may give you slightly different program versions (which in principle could lead to slightly different results).

Then activate the conda environment with:

    conda activate CoV_229E_antigenic_drift

and then run the pipeline in [Snakefile](Snakefile) with:

    snakemake --cores all

Here are the key steps performed by [Snakefile](Snakefile) and the key resulting outputs:

#### Get and parse Spike sequences
Download the sequences by accession and parse out Spikes from human isolates.
This is done by the Jupyter notebook [get_parse_spikes.ipynb](get_parse_spikes.ipynb), and the Markdown output of this notebook is in [results/get_parse_spikes.md](results/get_parse_spikes.md).
This step creates a FASTA file of **unaligned** Spike nucleotide sequences (in [results/spikes_unaligned_nt.fasta](results/spikes_unaligned_nt.fasta)) and a CSV with metadata for these spikes (in [results/spikes_metadata.csv](results/spikes_metadata.csv)).

#### Build codon and protein alignments of Spike
We use [mafft](https://mafft.cbrc.jp/alignment/software/) and the [HyPhy codon-aware alignment scripts](https://github.com/veg/hyphy-analyses/tree/master/codon-msa) to build codon and protein alignments of the Spikes.
This creates the alignment files [results/spikes_aligned_codon.fasta](results/spikes_aligned_codon.fasta) and [spikes_aligned_prots.fasta](spikes_aligned_prots.fasta).
This step is done by rules in [Snakefile](Snakefile).

#### Screen sequences for recombination

#### Infer phylogenetic tree
We infer a phylogenetic tree on the codon alignment using a codon-substitution model with [IQTREE](http://www.iqtree.org/).
The inferred tree is in [results/spikes_iqtree.treefile](results/spikes_iqtree.treefile).
This step is done by a rule in [Snakefile](Snakefile).
