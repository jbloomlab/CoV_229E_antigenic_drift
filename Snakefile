"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import os
import requests
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
hyphy_scripts = ['pre-msa.bf', 'post-msa.bf', 'lib/igscueal.bf']

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule all:
    """Target rule with final output files."""
    input:
        config['spikes_metadata'],
        os.path.join(config['results_dir'], 'get_parse_spikes.md'),
        config['spikes_unaligned_nt'],
        config['spikes_aligned_codon'],
        config['spikes_aligned_prot'],
        gard_recomb_json=config['gard_recomb_json'],
        gard_recomb_best=config['gard_recomb_best'],

rule gard_recomb_screen:
    """Use GARD to screen for recombination:
       https://academic.oup.com/bioinformatics/article/22/24/3096/208339
       https://link.springer.com/protocol/10.1007/978-1-4939-9074-0_14
    """
    input:
        spikes_aligned_codon=config['spikes_aligned_codon'],
    output:
        gard_recomb_json=config['gard_recomb_json'],
        gard_recomb_best=config['gard_recomb_best'],
    shell:
        """
        hyphy gard \
            --rv GDD \
            --rate-classes 3 \
            --alignment {input.spikes_aligned_codon} \
            --output {output.gard_recomb_json} \
            --output-lf {output.gard_recomb_best}
        """
        

rule build_codon_alignment:
    """Build codon alignment using `mafft` and `HyPhy` as here:
       https://github.com/veg/hyphy-analyses/tree/master/codon-msa
    """
    input:
        pre_msa=os.path.join(config['hyphy_scripts_dir'], 'pre-msa.bf'),
        post_msa=os.path.join(config['hyphy_scripts_dir'], 'post-msa.bf'),
        libscript=os.path.join(config['hyphy_scripts_dir'], 'lib/igscueal.bf'),
        spikes_unaligned_nt=config['spikes_unaligned_nt'],
    output:
        spikes_aligned_codon=config['spikes_aligned_codon'],
        spikes_aligned_prot=config['spikes_aligned_prot'],
        temp_prot=temp(config['spikes_unaligned_nt'] + '_protein.fas'),
        temp_nt=temp(config['spikes_unaligned_nt'] + '_nuc.fas'),
    shell:
        """
        hyphy {input.pre_msa} --input {input.spikes_unaligned_nt}
        mafft --auto {output.temp_prot} > {output.spikes_aligned_prot}
        hyphy {input.post_msa} \
            --protein-msa {output.spikes_aligned_prot} \
            --nucleotide-sequences {output.temp_nt} \
            --output {output.spikes_aligned_codon} \
            --compress No
        """

rule get_hyphy_codon_align_scripts:
    """Get scripts used by HyPhy to make codon alignments.
       https://github.com/veg/hyphy-analyses/tree/master/codon-msa
    """
    output:
        [os.path.join(config['hyphy_scripts_dir'], s) for s in
         ['pre-msa.bf', 'post-msa.bf', 'lib/igscueal.bf']]
    params:
        base_url='https://raw.githubusercontent.com/veg/hyphy-analyses/f298f6d8e29de2fd406c48747efe894b55d0c8ce/codon-msa/'
    run:
        for s in output:
            r = requests.get(os.path.join(
                            params.base_url,
                            os.path.relpath(s, config['hyphy_scripts_dir'])))
            with open(s, 'wb') as f:
                f.write(r.content)

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
