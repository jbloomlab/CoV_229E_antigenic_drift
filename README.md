# Antigenic drift of the human seasonal coronavirus 229E

Study by Laurel Kelnhofer-Millevolte, Rachel Eguia, Katharine Crawford, and Jesse Bloom

## Steps in phylogenetic analysis

### Obtaining Genbank accessions for sequence set
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
