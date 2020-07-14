# Get and parse Spike sequences
This Python Jupyter notebook identifies the human Spike sequences and writes them to a file.

## Imports and configuration
Import modules and read configuration file:


```python
import os
import re
import requests

import Bio.Entrez
import Bio.SeqIO

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import yaml
```

Read configuration:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)

# create directories
for dirname in [config['results_dir'], config['genbank_seqs_dir']]:
    os.makedirs(dirname, exist_ok=True)
```

## Download sequences
Read the CSV giving the Genbank accessions for all sequences, identify just the ones that are human sequences, then download all of them:


```python
# read the accession data frame
print(f"Reading accessions and metadata from {config['accessions']}\n")
acc_df = pd.read_csv(config['accessions'])
assert acc_df['Accession'].nunique() == len(acc_df)

# add additional manually specified metadata
print(f"Reading additional metadata from {config['addtl_accessions_metadata']}")
with open(config['addtl_accessions_metadata']) as f:
    addtl_metadata = (
            pd.DataFrame.from_dict(yaml.safe_load(f), orient='index')
            .rename_axis('Accession')
            .reset_index()
            )
assert set(addtl_metadata.columns).issubset(acc_df.columns)
assert set(addtl_metadata['Accession']).issubset(acc_df['Accession'])
# make sure additional metadata is only replacing NA values
replace_df = acc_df[addtl_metadata.columns].merge(addtl_metadata, on='Accession')
assert all(all(1 == replace_df[f"{col}_x"].isnull().astype(int) + replace_df[f"{col}_y"].isnull().astype(int))
           for col in addtl_metadata.columns if col != 'Accession')
# fill in additional metadata
addtl_acc = addtl_metadata['Accession'].tolist()
acc_df = (acc_df
          .set_index('Accession')
          .combine_first(addtl_metadata.set_index('Accession'))
          .reset_index()
          )

# identify and retain human sequences:
acc_df = acc_df.assign(is_human=lambda x: x['Host'] == 'Homo sapiens')
p = (ggplot(acc_df) +
     aes('Host', fill='is_human') +
     geom_bar() +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(0.3 * acc_df['Host'].nunique(), 2))
     )
_ = p.draw()
print(f"Retaining just the {len(acc_df.query('is_human == True'))} of {len(acc_df)} accessions that are human.")
acc_df = acc_df.query('is_human == True')

acc_to_file = {acc: os.path.join(config['genbank_seqs_dir'], f"{acc}.gb")
               for acc in acc_df['Accession']}
missing_accessions = {acc: f for acc, f in acc_to_file.items() if not os.path.isfile(f)}
print(f"There are {len(acc_to_file)} accessions in {config['accessions']}, of which "
      f"{len(missing_accessions)} are not already downloaded in {config['genbank_seqs_dir']}.")
if config['redownload_genbank']:
    to_download = acc_to_file
    print(f"Downloading all {len(to_download)} accessions (overwriting any existing ones).")
else:
    to_download = missing_accessions
    print(f"Downloading the {len(to_download)} accessions not already in {config['genbank_seqs_dir']}.")

Bio.Entrez.email = config['email']  # set e-mail address for downloads    
update_every = 25
for i, (acc, fname) in enumerate(to_download.items()):
    if i % update_every == 0:
        print(f"Progress: accession {i + 1} of {len(to_download)}")
    gb_text = Bio.Entrez.efetch(db='nucleotide',
                                id=acc,
                                rettype='gb',
                                retmode='text',
                                ).read()
    with open(fname, 'w') as f:
        f.write(gb_text)
```

    Reading accessions and metadata from data/NCBI_Virus_229E_accessions.csv
    
    Reading additional metadata from data/extra_229E_accessions_metadata.yaml
    Retaining just the 329 of 546 accessions that are human.
    There are 329 accessions in data/NCBI_Virus_229E_accessions.csv, of which 329 are not already downloaded in results/genbank_seqs.
    Downloading the 329 accessions not already in results/genbank_seqs.
    Progress: accession 1 of 329
    Progress: accession 26 of 329
    Progress: accession 51 of 329
    Progress: accession 76 of 329
    Progress: accession 101 of 329
    Progress: accession 126 of 329
    Progress: accession 151 of 329
    Progress: accession 176 of 329
    Progress: accession 201 of 329
    Progress: accession 226 of 329
    Progress: accession 251 of 329
    Progress: accession 276 of 329
    Progress: accession 301 of 329
    Progress: accession 326 of 329



![png](get_parse_spikes_files/get_parse_spikes_6_1.png)


## Parse Spike coding sequences
Make a data frame that has a row for each CDS in each Genbank accession (these accessions often represent many CDSs), identify what products these CDSs correspond to, and extract the CDS sequence.
Furthermore, use the CDS product name to determine if it is Spike.
After doing this, we plot the number of CDSs for each product name, colored by whether they are Spike.
Look at this plot to make sure that everything colored as Spike indeed appears to be Spike, and everything **not** colored as Spike is not Spike:


```python
cds_df = (
    acc_df
    .assign(SeqRecord=lambda x: (x['Accession']
                                 .map(acc_to_file)
                                 .apply(Bio.SeqIO.read, format='gb')
                                 ),
            cds_features=lambda x: x['SeqRecord'].map(lambda s: [f for f in s.features if f.type == 'CDS']),
            n_cds_features=lambda x: x['cds_features'].map(len)
            )
    .query('n_cds_features > 0')
    .explode('cds_features')
    .rename(columns={'cds_features': 'cds_feature'})
    .assign(product=lambda x: x['cds_feature'].map(lambda f: f.qualifiers['product'][0]),
            sequence=lambda x: x.apply(lambda r: r['cds_feature'].extract(r['SeqRecord']).seq,
                                       axis=1),
            sequence_length=lambda x: x['sequence'].map(len),
            is_spike=lambda x: x['product'].str.match('^(S|s)'),
            )
    )

p = (ggplot(cds_df) +
     aes('product', fill='is_spike') +
     geom_bar() +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(0.2 * cds_df['product'].nunique(), 2.5))
     )
_ = p.draw()
```


![png](get_parse_spikes_files/get_parse_spikes_8_0.png)


Now get a data frame that is just Spike, filtering just to those of valid lengths, and checking for problematic features like premature stop codons or ambiguous nucleotides:


```python
spike_df = (
    cds_df
    .query('is_spike == True')
    .assign(protein=lambda x: x['sequence'].map(lambda s: s[: (len(s) // 3) * 3].translate()).map(str),
            protein_length=lambda x: x['protein'].map(len),
            valid_length=lambda x: x['protein_length'].map(lambda n: config['prot_length_range'][0] <= n
                                                                     <= config['prot_length_range'][1]),
            premature_stop=lambda x: x['protein'].str[: -1].str.contains('*', regex=False),
            sequence=lambda x: x['sequence'].map(str),
            ambiguous_nts=lambda x: x['sequence'].str.contains('[^ACGT]'),
            country=lambda x: x['Geo_Location'].str.split(':').map(lambda y: y[0]),
            collection_date=lambda x: pd.to_datetime(x['Collection_Date']),
            no_collection_date=lambda x: x['collection_date'].isnull(),
            )
    )

print('Here is plot of number of Spikes at each protein length; only keeping those with valid length:')
p = (ggplot(spike_df.assign(protein_length=lambda x: pd.Categorical(x['protein_length']))) +
     aes('protein_length', fill='valid_length') +
     geom_bar() +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(0.3 * spike_df['protein_length'].nunique(), 2.5),
           )
     ).draw()
display(p)
plt.close(p)

spike_df = spike_df.query('valid_length')

# columns we keep to print
cols_of_interest = ['Accession', 'Authors', 'Geo_Location', 'country', 'Isolation_Source',
                    'collection_date', 'GenBank_Title', 'sequence', 'protein']

for filter_criteria in ['premature_stop', 'ambiguous_nts', 'no_collection_date']:
    df_fail = spike_df.query(filter_criteria)
    print(f"\nRemoving {len(df_fail)} of {len(spike_df)} Spikes for failing {filter_criteria} filter.")
    if len(df_fail):
        print('Here are the sequences being filtered:')
        display(df_fail[cols_of_interest])
        spike_df = spike_df.query(f"not {filter_criteria}")
        
# add strain names
def get_strain_name(genbank_title):
    """Parse meaningful strain name from `GenBank_Title`"""
    m = re.match('Human coronavirus 229E (?:isolate|strain|S gene for spike protein, complete cds, isolate\:) '
                 '(?:HCoV_229E/|HCoV\-229E/|229E/human/|229E/)?'
                 '([^\s,]+)'
                 '(?: spike|, |$)',
                 genbank_title)
    if m is None:
        raise ValueError(f"cannot match GenBank_Title:\n{genbank_title}")
    return m.group(1)

# add nice names and sort by date, and get columns of interest
cols_of_interest.insert(0, 'strain_name')
cols_of_interest.append('year')
spike_df = (
    spike_df
    .assign(strain_name=lambda x: x['GenBank_Title'].map(get_strain_name),
            year=lambda x: x['collection_date'].dt.year,
            )
    .sort_values('collection_date')
    .reset_index()
    [cols_of_interest]
    )
assert len(spike_df) == spike_df['strain_name'].nunique()

print(f"\nOverall, retained {len(spike_df)} Spikes.")
```

    Here is plot of number of Spikes at each protein length; only keeping those with valid length:



![png](get_parse_spikes_files/get_parse_spikes_10_1.png)


    
    Removing 0 of 87 Spikes for failing premature_stop filter.
    
    Removing 2 of 87 Spikes for failing ambiguous_nts filter.
    Here are the sequences being filtered:



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Accession</th>
      <th>Authors</th>
      <th>Geo_Location</th>
      <th>country</th>
      <th>Isolation_Source</th>
      <th>collection_date</th>
      <th>GenBank_Title</th>
      <th>sequence</th>
      <th>protein</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>167</th>
      <td>JX503060</td>
      <td>Farsani,S.M., Dijkman,R., Jebbink,M.F., Goosse...</td>
      <td>Netherlands</td>
      <td>Netherlands</td>
      <td>oronasopharynx</td>
      <td>2010-01-01</td>
      <td>Human coronavirus 229E isolate 0349, complete ...</td>
      <td>ATGTTTGTTTTACTTGTTGCATATGCCTTGTTGCATATTGCTGGTT...</td>
      <td>MFVLLVAYALLHIAGCQTTNGTNTSHSVCNGCVGHSENVFAVESGG...</td>
    </tr>
    <tr>
      <th>533</th>
      <td>MT438699</td>
      <td>Dinwiddie,D.L., Dehority,W.N., Schwalm,K.C., K...</td>
      <td>USA: Little Rock, Arkansas</td>
      <td>USA</td>
      <td>oronasopharynx</td>
      <td>2017-02-27</td>
      <td>Human coronavirus 229E isolate HCoV-229E/USA/A...</td>
      <td>ATGTTTGTTTTACTTGTTGCATATGCCTTGTTGCATATTGCTGGTT...</td>
      <td>MFVLLVAYALLHIAGCQTTNGTNTSHSVCNGCVGHSENVFAVESGG...</td>
    </tr>
  </tbody>
</table>
</div>


    
    Removing 3 of 85 Spikes for failing no_collection_date filter.
    Here are the sequences being filtered:



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Accession</th>
      <th>Authors</th>
      <th>Geo_Location</th>
      <th>country</th>
      <th>Isolation_Source</th>
      <th>collection_date</th>
      <th>GenBank_Title</th>
      <th>sequence</th>
      <th>protein</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>189</th>
      <td>KF293664</td>
      <td>Lundin,A., Dijkman,R., Bergstrom,T., Kann,N., ...</td>
      <td>Sweden</td>
      <td>Sweden</td>
      <td>NaN</td>
      <td>NaT</td>
      <td>Human coronavirus 229E clone p0, partial genome</td>
      <td>ATGTTTGTTTTGCTTGTTGCATATGCCTTGTTGCATATTGCTGGTT...</td>
      <td>MFVLLVAYALLHIAGCQTTNGLNTSYSVCNGCVGYSENVFAVESGG...</td>
    </tr>
    <tr>
      <th>190</th>
      <td>KF293665</td>
      <td>Lundin,A., Dijkman,R., Bergstrom,T., Kann,N., ...</td>
      <td>Sweden</td>
      <td>Sweden</td>
      <td>NaN</td>
      <td>NaT</td>
      <td>Human coronavirus 229E clone mock-p11, partial...</td>
      <td>ATGTTTGTTTTGCTTGTTGCATATGCCTTGTTGCATATTGCTGGTT...</td>
      <td>MFVLLVAYALLHIAGCQTTNGLNTSYSVCNGCVGYSENVFAVESGG...</td>
    </tr>
    <tr>
      <th>191</th>
      <td>KF293666</td>
      <td>Lundin,A., Dijkman,R., Bergstrom,T., Kann,N., ...</td>
      <td>Sweden</td>
      <td>Sweden</td>
      <td>NaN</td>
      <td>NaT</td>
      <td>Human coronavirus 229E clone K22-p11, partial ...</td>
      <td>ATGTTTGTTTTGCTTGTTGCATATGCCTTGTTGCATATTGCTGGTT...</td>
      <td>MFVLLVAYALLHIAGCQTTNGLNTSYSVCNGCVGYSENVFAVESGG...</td>
    </tr>
  </tbody>
</table>
</div>


    
    Overall, retained 82 Spikes.


## Examine distribution of sequences by year and country
Plot Spike sequence isolation year and location (country):


```python
p = (ggplot(spike_df) +
     aes('year') +
     geom_histogram(binwidth=1) +
     facet_wrap('~ country', nrow=1) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2 * spike_df['country'].nunique(), 1.5),
           )
     )

_ = p.draw()
```


![png](get_parse_spikes_files/get_parse_spikes_12_0.png)


## Write Spike sequences and metadata to files
Write the Spike sequences to a FASTA file with the strain names as headers, and write the metadata for each strain to a CSV file:


```python
print(f"Writing the Spike nucleotide sequences to {config['spikes_unaligned_nt']} and the "
      f"metadata to {config['spikes_metadata']}")

(spike_df
 .drop(columns=['sequence', 'protein'])
 .to_csv(config['spikes_metadata'], index=False)
 )

with open(config['spikes_unaligned_nt'], 'w') as f:
    for row in spike_df.itertuples():
        f.write(f">{row.strain_name}\n{row.sequence}\n")
```

    Writing the Spike nucleotide sequences to results/spikes_unaligned_nt.fasta and the metadata to results/spikes_metadata.csv
