"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import os
import subprocess

# Configuration ---------------------------------------------------------------
configfile: 'config.yaml'

# Functions -------------------------------------------------------------------

def run_nb_to_md(nb):
    """Run Jupyter notebook and put markdown results in `results_dir`."""
    subprocess.check_call(['jupyter', 'nbconvert',
                           '--to', 'notebook',
                           '--execute',
                           '--inplace',
                           '--ExecutePreprocessor.timeout=-1',
                           nb,
                           ])
    subprocess.check_call(['jupyter', 'nbconvert',
                           '--output-dir', config['results_dir'],
                           '--to', 'markdown',
                           nb,
                           ])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule all:
    """Target rule with final output files."""
    input:
        config['spikes_metadata'],
        os.path.join(config['results_dir'], 'get_parse_spikes.md'),

rule get_parse_spikes:
    """Get and parse human Spike sequences from Genbank."""
    input:
        config['accessions'],
        config['addtl_accessions_metadata'],
    output:
        config['spikes_unaligned_nt'],
        config['spikes_metadata'],
        os.path.join(config['results_dir'], 'get_parse_spikes.md'),
    run:
        run_nb_to_md('get_parse_spikes.ipynb')
