# Antigenic drift of the human seasonal coronavirus 229E

Study by Laurel Kelnhofer-Millevolte, Rachel Eguia, Katharine Crawford, and Jesse Bloom

# Steps in phylogenetic analysis

1. Manually download 229E sequences from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/).
   To do this, go to the [NCBI Virus website](https://www.ncbi.nlm.nih.gov/labs/virus/), click *Search by virus*, and then enter *Human coronavirus 229E* (taxid:1137) in the *Search by virus name or taxonomy* box and hit return.
   Then select all the sequences, click *Download*, select *Coding region* under *Sequence data (FASTA format)*, hit *Next*, click *Download Selected Records* (they should all be selected), hit *Next* again and set up the FASTA definition line to include:
     - CDS id
     - GenBank Titel
     - Host
     - Species
     - Geo Location

   and then save the file as [data/Genbank_229E_sequences.fasta](data/Genbank_229E_sequences.fasta).
   The current file was downloaded on July 10, 2020.

