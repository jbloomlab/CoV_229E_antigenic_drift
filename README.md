# Antigenic drift of the human seasonal coronavirus 229E

Study by Laurel Kelnhofer-Millevolte, Rachel Eguia, Katharine Crawford, and Jesse Bloom

The analysis is divided into some manual steps that provide the input data, followed by an automated computational analysis to prodcue the rest of the results from input data.

## Manual steps
These manual steps are used to generate the input data in the [./data/](data) subdirectory.

### Genbank accessions and extra metadata
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

In addition, we want to specify certain accessions to ensure they are included, excluded, or specially annotated (for instance because they are being used in experiments).
The file [data/accessions_to_include_exclude_annotate.yaml](data/accessions_to_include_exclude_annotate.yaml) specifies this information; the file should be self-explanatory.

### Domain-annotated Spike
For domain-level analysis, we use a GenePept (`*.gp`) file that contains the Spike protein that was used in the PDB structure reported by [Li et al (2019)](https://elifesciences.org/articles/51230) (this is PDB [6u7h](https://www.rcsb.org/structure/6U7H)).
The Genbank protein accession for the Spike used in that study is [AAK32191](https://www.ncbi.nlm.nih.gov/protein/AAK32191).
This accession was downloaded in GenePept format, and then manually annotated to describe key domains / motifs (using the definitions in [Li et al (2019)](https://elifesciences.org/articles/51230)) to create [data/AAK32191_hand_annotated.gp](data/AAK32191_hand_annotated.gp).

## Steps run automatically by Snakemake pipeline
The remaining steps are run automatically using [snakemake](https://snakemake.readthedocs.io/) to run [Snakefile](Snakefile), which reads its configuration from [config.yaml](config.yaml).
The results of the automated steps are placed in [./results/](results).

To run these steps, first build the conda environment, which installs the necessary programs.
First [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
Then build the environment.
You can probably install the full version-pinned environment in [environment_pinned.yml](environment_pinned.yml) with:

    conda env create -f environment_pinned.yml

If this fails (which is possible if you are using a different operating system), then install the less-version-unpinned environment in [environment_unpinned.yml](environment_unpinned.yml), which may give you slightly different program versions (which in principle could lead to slightly different results).

Then activate the conda environment with:

    conda activate CoV_229E_antigenic_drift

and then run the pipeline in [Snakefile](Snakefile) with:

    snakemake --cores all

Here are the key steps performed by the automated analysis in [Snakefile](Snakefile) and the key outputs:

### Get and parse Spike sequences
Download the sequences by accession and parse out Spikes from human isolates.
This is done by the Jupyter notebook [get_parse_spikes.ipynb](get_parse_spikes.ipynb), and the Markdown output of this notebook is in [results/get_parse_spikes.md](results/get_parse_spikes.md).
This step creates a CSV with metadata for the spikes in [results/spikes_metadata.csv](results/spikes_metadata.csv).

### Build codon and protein alignments of Spike
We use [mafft](https://mafft.cbrc.jp/alignment/software/) and [a custom Python script](prot_to_codon_alignment.py) to build protein and codon alilgnments of the Spikes.
This creates the alignment files [results/spikes_aligned_codon.fasta](results/spikes_aligned_codon.fasta) and [results/spikes_aligned_prots.fasta](results/spikes_aligned_prots.fasta).

### Screen sequences for recombination

### Infer phylogenetic tree and root and time-scale it.
We infer a phylogenetic tree from the codon alignment using a codon-substitution model with [IQTREE](http://www.iqtree.org/).
We then use [treetime](https://treetime.readthedocs.io/) to root the resulting tree, and also create a time-scaled tree in which the branch lengths are adjusted to be in units of time.
The time-scaling assumes that branch lengths are proportional to time, so you should check this assumption by looking at the root-to-tip regression in [results/timetree/root_to_tip_regression.pdf](results/timetree/root_to_tip_regression.pdf).
The tree with branch lengths still in units of codon substitutions per site is in [results/timetree/divergence_tree.newick](results/timetree/divergence_tree.newick), and the tree with branch lengths scaled by time is in [results/timetree/timetree.newick](results/timetree/timetree.newick).

### Draw phylogenetic trees
We use [ete3](http://etetoolkit.org/) to make nice drawings of the phylogenetic tree (both the divergence and time tree).
This is done by the Jupyter notebook [draw_trees.ipynb](draw_trees.ipynb), and the Markdown output of this notebook is in [results/draw_trees.md](results/draw_trees.md).
The time-scaled tree drawing in [results/tree/timetree.pdf](results/tree/timetree.pdf) and the divergence-scaled tree drawing is in [results/tree/divergence_tree.pdf](results/tree/divergence_tree.pdf).
A legend that maps the tip coloring to country of isolation is in [results/tree/legend.pdf](results/tree/legend.pdf).

### Variability of Spike domains
We compute the variability of each site in Spike over the sequence alignment, and plot this on the domain structure (using the residue number and domain definitions defined in [data/AAK32191_hand_annotated.gp](data/AAK32191_hand_annotated.gp)).
This is done by the Jupyter notebook [analyze_variation.ipynb](analyze_variation.ipynb), and the Markdown output of this notebook is in [results/analyze_variation.md](results/analyze_variation.md).
The output is the schematic [results/variation_analysis/domain_schematic.pdf](results/variation_analysis/domain_schematic.pdf).

In addition, the created file [results/variation_analysis/site_variability.csv](results/variation_analysis/site_variability.csv) gives the variability at each site.
This file along with the PDB at [results/pdbs/6u7h.pdb](results/pdbs/6u7h.pdb) and the metadata file at [results/variation_analysis/dms_view_metadata.md](results/variation_analysis/dms_view_metadata.md) can be used to visualize the variability on the protein structure.

### Design Spikes for experiments
We design the sequences of the Spikes that will be used for the experiments (they have a C-terminal tail deletion to improve lentiviral pseudotyping).
This design is done by the Juyter notebook [seqs_for_expts.ipynb](seqs_for_expts.ipynb), and the Markdown output of this notebook is in [results/seqs_for_expts.md](results/seqs_for_expts.md).
The amino-acid Levenshtein distance (essentially number of differences including gaps) among these sequences is plotted in [results/seqs_for_expts/n_aa_diffs.pdf](results/seqs_for_expts/n_aa_diffs.pdf)
