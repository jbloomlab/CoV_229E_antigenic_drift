"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import os
import requests
import subprocess

import Bio.Phylo

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
        config['spikes_unaligned_nt'],
        config['spikes_unaligned_prot'],
        config['spikes_aligned_prot'],
        config['spikes_aligned_codon'],
        config['divergencetree'],
        config['timetree'],
        config['root_to_tip'],
        os.path.join(config['results_dir'], 'draw_trees.md'),
        config['tree_legend'],
        config['divergencetree_image'],
        config['timetree_image'],
        os.path.join(config['pdb_dir'], f"{config['pdb_id']}.pdb"),
        config['site_variation'],
        config['dms_view_metadata'],
        config['variation_domain_schematic'],
        os.path.join(config['results_dir'], 'analyze_variation.md'),
#        gard_recomb_json=config['gard_recomb_json'],
#        gard_recomb_best=config['gard_recomb_best'],

rule variability_analysis:
    input:
        os.path.join(config['pdb_dir'], f"{config['pdb_id']}.pdb"),
        config['domain_annotated_spike'],
        config['spikes_unaligned_prot']
    output:
        config['site_variation'],
        config['dms_view_metadata'],
        config['variation_domain_schematic'],
        os.path.join(config['results_dir'], 'analyze_variation.md'),
    run:
        run_nb_to_md('analyze_variation.ipynb')
        
rule get_pdb:
    output:
        pdbfile=os.path.join(config['pdb_dir'], "{pdb_id}.pdb")
    run:
        pdb_id = os.path.splitext(os.path.basename(output.pdbfile))[0]
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        r = requests.get(url)
        with open(output.pdbfile, 'wb') as f:
            f.write(r.content)

rule draw_tree:
    """Draw the phylogenetic trees."""
    input:
        config['accessions_special'],
        config['spikes_metadata'],
        config['divergencetree'],
        config['timetree'],
    output:
        config['tree_legend'],
        config['divergencetree_image'],
        config['timetree_image'],
        os.path.join(config['results_dir'], 'draw_trees.md'),
    run:
        run_nb_to_md('draw_trees.ipynb')

rule nexus_to_newick:
    """Convert tree from Nexus to Newick format."""
    input:
        nexus="{treename}.nexus"
    output:
        newick="{treename}.newick"
    run:
        Bio.Phylo.convert(input.nexus, 'nexus', output.newick, 'newick')

rule timetree:
    """Run ``treetime`` to root and make time-scaled tree."""
    input:
        aln=config['spikes_aligned_codon'],
        dates=config['spikes_metadata'],
        tree='results/iqtree/spikes.treefile',
    output:
        os.path.splitext(config['divergencetree'])[0] + '.nexus',
        os.path.splitext(config['timetree'])[0] + '.nexus',
        config['root_to_tip'],
    params:
        outdir=os.path.dirname(config['timetree'])
    shell:
        """
        treetime \
            --aln {input.aln} \
            --dates {input.dates} \
            --tree {input.tree} \
            --outdir {params.outdir} \
            --verbose 1
        """

rule build_iqtree:
    """Use IQTREE to infer a phylogenetic tree."""
    input:
        spikes_aligned_codon=config['spikes_aligned_codon'],
    output:
        spikes_iqtree='results/iqtree/spikes.treefile',
    params:
        prefix='results/iqtree/spikes',
    shell:
        """
        iqtree \
            -s {input.spikes_aligned_codon} \
            -st CODON \
            -m MGK+G+F3X4 \
            -pre {params.prefix}
        """

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
    """Create codon alignment from protein alignment and nucleotide sequences."""
    input:
        aligned_prot=config['spikes_aligned_prot'],
        nt_seqs=config['spikes_unaligned_nt'],
    output:
        aligned_codon=config['spikes_aligned_codon']
    shell:
        """
        python prot_to_codon_alignment.py \
            --prot_aln {input.aligned_prot} \
            --nt_seqs {input.nt_seqs} \
            --codon_aln {output.aligned_codon}
        """

rule align_spike_prots:
    """Use ``mafft`` to create a protein alignment."""
    input:
        unaligned_prot=config['spikes_unaligned_prot'],
    output:
        aligned_prot=config['spikes_aligned_prot'],
    shell:
        "mafft --auto {input.unaligned_prot} > {output.aligned_prot}"

rule get_parse_spikes:
    """Get and parse human Spike sequences from Genbank."""
    input:
        config['accessions'],
        config['addtl_accessions_metadata'],
        config['accessions_special'],
    output:
        config['spikes_unaligned_nt'],
        config['spikes_unaligned_prot'],
        config['spikes_metadata'],
        os.path.join(config['results_dir'], 'get_parse_spikes.md'),
    run:
        run_nb_to_md('get_parse_spikes.ipynb')
