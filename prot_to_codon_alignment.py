"""Creates codon alignment from protein alignment and nucleotide sequences.

For usage, type::

    python prot_to_codon_alignment.py --help

"""


import argparse
import os

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord


def main():

    parser = argparse.ArgumentParser(
                description=('Codon alignment from protein alignment and '
                             'nucleotide sequences.')
                )
    parser.add_argument('--prot_aln',
                        required=True,
                        help='Input protein alignment (FASTA).'),
    parser.add_argument('--nt_seqs',
                        required=True,
                        help='Input nucleotide sequences (FASTA).'),
    parser.add_argument('--codon_aln',
                        required=True,
                        help='Name of created codon alignment (FASTA).'),
    args = parser.parse_args()

    prot_aln = list(Bio.SeqIO.parse(args.prot_aln, 'fasta'))
    print(f"Read {len(prot_aln)} aligned proteins from {args.prot_aln}")
    prot_len = len(prot_aln[0])
    if any(prot_len != len(p) for p in prot_aln):
        raise ValueError('the aligned proteins not all of same length')
    print(f"The aligned proteins are of length {prot_len} residues.")

    nt_seqs = list(Bio.SeqIO.parse(args.nt_seqs, 'fasta'))
    print(f"Read {len(nt_seqs)} nucleotide sequences from {args.nt_seqs}")

    # some error checks
    if len(prot_aln) != len(nt_seqs):
        raise ValueError('inconsistent number of sequences in protein '
                         'alignment and nucleotide sequences')
    for i, (p, nt) in enumerate(zip(prot_aln, nt_seqs)):
        if p.id != nt.id or p.description != nt.description:
            raise ValueError(f"entry {i + 1} differs in protein alignment and "
                             f"nucleotide sequences:\n{p.id} {p.description}"
                             f"\n{nt.id} {nt.description}")
        if len(nt) != 3 * (len(p) - p.seq.count('-')):
            raise ValueError(f"inconsistent length: {p.id} {p.description}")

    # now get aligned codon sequences
    codon_aln = []
    for p, nt in zip(prot_aln, nt_seqs):
        c = []
        r = 0
        ntseq = str(nt.seq)
        for aa in p.seq:
            if aa != '-':
                c.append(ntseq[3 * r: 3 * r + 3])
                r += 1
            else:
                c.append('---')
        c = ''.join(c)
        assert len(c) - c.count('-') == len(nt) == 3 * (len(p) - p.seq.count('-'))
        codon_aln.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(c),
                                                 id=nt.id,
                                                 name=nt.name,
                                                 description=nt.description))
    assert len(codon_aln) == len(nt_seqs) == len(prot_aln)

    print(f"Writing codon alignment to {args.codon_aln}")
    if os.path.dirname(args.codon_aln):
        os.makedirs(os.path.dirname(args.codon_aln), exist_ok=True)
    Bio.SeqIO.write(codon_aln, args.codon_aln, 'fasta')


if __name__ == '__main__':
    main()
