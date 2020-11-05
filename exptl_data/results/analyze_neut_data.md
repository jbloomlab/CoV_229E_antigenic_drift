# Analyze neutralization data
This Python Jupyter notebook analyzes the neutralization data.

Import Python modules.
We use [neutcurve](https://jbloomlab.github.io/neutcurve/) to plot the neutralization curves:


```python
import re
import warnings

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import neutcurve
from neutcurve.colorschemes import CBPALETTE
from neutcurve.colorschemes import CBMARKERS

import pandas as pd

from plotnine import *

print(f"Using `neutcurve` version {neutcurve.__version__}")
```

    Using `neutcurve` version 0.5.0


Specify input / output files:


```python
# input files
fracinfect_file = 'results/fracinfect.csv'
serum_info_file = 'serum_info.csv'

# output files
all_replicate_curves = 'results/all_neut_replicates.pdf'
all_neut_by_sera_curves = 'results/all_neut_by_sera.pdf'
all_fit_params = 'results/all_fit_params.csv'
all_neut_titers = 'results/all_neut_titers.csv'
```

Read in the neutralization data, dropping sera labeled as "Nothing":


```python
print(f"Reading neutralization data from {fracinfect_file}")
fracinfect = pd.read_csv(fracinfect_file)

# Input data restarts replicate numbers for each new date. Instead label replicates as:
#  - replicate_on_date: number of the replicate on that specific date
#  - replicate_with_date: number of replicate on that date suffixed by date
#  - replicate_all_dates: number replicates sequentially across all dates
fracinfect = (
    fracinfect
    .query('serum != "Nothing"')
    .assign(replicate_with_date=lambda x: x['replicate'].astype(str) +
                                          ' (' + x['date'] + ')')
    .rename(columns={'replicate': 'replicate_on_date'})
    )
fracinfect = (
    fracinfect
    .merge(fracinfect
           .sort_values('date')
           [['serum', 'virus', 'replicate_with_date']]
           .drop_duplicates()
           .assign(replicate_all_dates=lambda x: x.groupby(['serum', 'virus'])
                                                  ['replicate_with_date']
                                                  .transform('cumcount') + 1
                   ),
            how='left', on=['serum', 'virus', 'replicate_with_date'], validate='many_to_one',
            )
    )

# make sure unique reading for each virus / serum / replicate / date
assert len(fracinfect) == len(fracinfect.groupby(['serum',
                                                  'virus',
                                                  'replicate_all_dates',
                                                  'concentration',
                                                  ]))

# show first few lines of data frame
display(HTML(fracinfect.head().to_html(index=False)))
```

    Reading neutralization data from results/fracinfect.csv



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>replicate_on_date</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>date</th>
      <th>replicate_with_date</th>
      <th>replicate_all_dates</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>FH007TR</td>
      <td>229E-1992</td>
      <td>1</td>
      <td>0.100000</td>
      <td>0.4530</td>
      <td>2020-10-22</td>
      <td>1 (2020-10-22)</td>
      <td>1</td>
    </tr>
    <tr>
      <td>FH007TR</td>
      <td>229E-1992</td>
      <td>1</td>
      <td>0.033330</td>
      <td>0.7756</td>
      <td>2020-10-22</td>
      <td>1 (2020-10-22)</td>
      <td>1</td>
    </tr>
    <tr>
      <td>FH007TR</td>
      <td>229E-1992</td>
      <td>1</td>
      <td>0.011110</td>
      <td>0.8930</td>
      <td>2020-10-22</td>
      <td>1 (2020-10-22)</td>
      <td>1</td>
    </tr>
    <tr>
      <td>FH007TR</td>
      <td>229E-1992</td>
      <td>1</td>
      <td>0.003704</td>
      <td>0.9029</td>
      <td>2020-10-22</td>
      <td>1 (2020-10-22)</td>
      <td>1</td>
    </tr>
    <tr>
      <td>FH007TR</td>
      <td>229E-1992</td>
      <td>1</td>
      <td>0.001235</td>
      <td>1.2490</td>
      <td>2020-10-22</td>
      <td>1 (2020-10-22)</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Get information on sera:


```python
serum_info = (
    pd.read_csv(serum_info_file)
    .assign(collection_date=lambda x: pd.to_datetime(x['collection_date']),
            collection_year=lambda x: x['collection_date'].dt.year +
                                      (x['collection_date'].dt.dayofyear - 1) / 365,
            )
    )
    
assert len(serum_info) == serum_info['serum'].nunique()
print(f"Read information for {len(serum_info)} sera")

sera_lacking_info = set(fracinfect['serum']) - set(serum_info['serum'])
if sera_lacking_info:
    raise ValueError(f"lacking information for these sera: {sera_lacking_info}")

# show first few lines of sera information data frame:
display(HTML(serum_info.head().to_html(index=False)))
```

    Read information for 28 sera



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>collection_date</th>
      <th>age</th>
      <th>category</th>
      <th>notes</th>
      <th>collection_year</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>SD85_1</td>
      <td>1985-01-03</td>
      <td>need from Terry</td>
      <td>human</td>
      <td>sera</td>
      <td>1985.005479</td>
    </tr>
    <tr>
      <td>SD85_2</td>
      <td>1985-01-08</td>
      <td>need from Terry</td>
      <td>human</td>
      <td>sera</td>
      <td>1985.019178</td>
    </tr>
    <tr>
      <td>SD85_3</td>
      <td>1985-04-10</td>
      <td>need from Terry</td>
      <td>human</td>
      <td>sera</td>
      <td>1985.271233</td>
    </tr>
    <tr>
      <td>SD85_4</td>
      <td>1985-04-11</td>
      <td>need from Terry</td>
      <td>human</td>
      <td>sera</td>
      <td>1985.273973</td>
    </tr>
    <tr>
      <td>SD85_5</td>
      <td>1985-04-18</td>
      <td>need from Terry</td>
      <td>human</td>
      <td>sera</td>
      <td>1985.293151</td>
    </tr>
  </tbody>
</table>


Use [neutcurve](https://jbloomlab.github.io/neutcurve/) to fit neutralization curves to all of the data:


```python
fits = neutcurve.curvefits.CurveFits(
            data=fracinfect,
            replicate_col='replicate_all_dates',
            fixbottom=0,
            fixtop=1,
            )
```

Plot all curves for all replicates of all virus / serum combinations:


```python
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # ignore fitting warnings
    fig, _ = fits.plotReplicates(ncol=8,
                                 legendtitle='replicate',
                                 xlabel='serum dilution',
                                 )
    
print(f"Saving plot to {all_replicate_curves}\n")
fig.savefig(all_replicate_curves)
fig.tight_layout()
display(fig)
plt.close(fig)
```

    Saving plot to results/all_neut_replicates.pdf
    



    
![png](analyze_neut_data_files/analyze_neut_data_11_1.png)
    


Make a plot showing all viruses against each sera:


```python
fig, _ = fits.plotSera(xlabel='serum dilution',
                       ncol=6,
                       legendtitle='virus')

print(f"Saving plot to {all_neut_by_sera_curves}\n")
fig.savefig(all_neut_by_sera_curves)
fig.tight_layout()
display(fig)
plt.close(fig)
```

    /fh/fast/bloom_j/software/miniconda3/envs/CoV_229E_antigenic_drift/lib/python3.8/site-packages/neutcurve/hillcurve.py:689: RuntimeWarning: invalid value encountered in power


    Saving plot to results/all_neut_by_sera.pdf
    



    
![png](analyze_neut_data_files/analyze_neut_data_13_2.png)
    


Write all of the fit parameters to a file:


```python
print(f"Writing all fit parameters to {all_fit_params}; first few lines also printed below:")

display(HTML(fits.fitParams().head().to_html(index=False)))

fits.fitParams().to_csv(all_fit_params, index=False)
```

    Writing all fit parameters to results/all_fit_params.csv; first few lines also printed below:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>FH007TR</td>
      <td>229E-1992</td>
      <td>average</td>
      <td>2</td>
      <td>0.100000</td>
      <td>lower</td>
      <td>&gt;0.1</td>
      <td>0.103966</td>
      <td>1.414215</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <td>FH007TR</td>
      <td>229E-1984</td>
      <td>average</td>
      <td>2</td>
      <td>0.008009</td>
      <td>interpolated</td>
      <td>0.00801</td>
      <td>0.008009</td>
      <td>1.074526</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <td>FH007TR</td>
      <td>229E-2001</td>
      <td>average</td>
      <td>2</td>
      <td>0.100000</td>
      <td>lower</td>
      <td>&gt;0.1</td>
      <td>0.102537</td>
      <td>1.004779</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <td>FH008WC</td>
      <td>229E-1992</td>
      <td>average</td>
      <td>2</td>
      <td>0.084957</td>
      <td>interpolated</td>
      <td>0.085</td>
      <td>0.084957</td>
      <td>1.268635</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <td>FH008WC</td>
      <td>229E-1984</td>
      <td>average</td>
      <td>2</td>
      <td>0.001461</td>
      <td>interpolated</td>
      <td>0.00146</td>
      <td>0.001461</td>
      <td>0.765615</td>
      <td>1</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Get neutralization titers, and merge them with information on the sera:


```python
neut_titers = (
    fits.fitParams()
    .assign(neut_titer=lambda x: 1 / x['ic50'])
    .assign(is_upper_bound=lambda x: x['ic50_bound'].map({'lower': True,
                                                          'interpolated': False}))
    [['serum', 'virus', 'neut_titer', 'is_upper_bound']]
    .merge(serum_info[['serum', 'collection_date', 'collection_year', 'age']],
           how='left', validate='many_to_one', on='serum')
    )

print(f"Writing neut titers to {all_neut_titers}; first few lines also printed below:")

display(HTML(neut_titers.head().to_html(index=False, float_format='%.1f')))

neut_titers.to_csv(all_neut_titers, index=False, float_format='%.1f')
```

    Writing neut titers to results/all_neut_titers.csv; first few lines also printed below:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>neut_titer</th>
      <th>is_upper_bound</th>
      <th>collection_date</th>
      <th>collection_year</th>
      <th>age</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>FH007TR</td>
      <td>229E-1992</td>
      <td>10.0</td>
      <td>True</td>
      <td>1989-08-29</td>
      <td>1989.7</td>
      <td>28</td>
    </tr>
    <tr>
      <td>FH007TR</td>
      <td>229E-1984</td>
      <td>124.9</td>
      <td>False</td>
      <td>1989-08-29</td>
      <td>1989.7</td>
      <td>28</td>
    </tr>
    <tr>
      <td>FH007TR</td>
      <td>229E-2001</td>
      <td>10.0</td>
      <td>True</td>
      <td>1989-08-29</td>
      <td>1989.7</td>
      <td>28</td>
    </tr>
    <tr>
      <td>FH008WC</td>
      <td>229E-1992</td>
      <td>11.8</td>
      <td>False</td>
      <td>1989-02-06</td>
      <td>1989.1</td>
      <td>unknown</td>
    </tr>
    <tr>
      <td>FH008WC</td>
      <td>229E-1984</td>
      <td>684.4</td>
      <td>False</td>
      <td>1989-02-06</td>
      <td>1989.1</td>
      <td>unknown</td>
    </tr>
  </tbody>
</table>


Now for plotting annotate neut titers with information about:
 - virus year
 - titer of serum against the 1984 virus


```python
def get_virus_year(virus):
    m = re.fullmatch('229E\-(?P<year>\d{4})', virus)
    if not m:
        raise ValueError(f"cannot match virus year in {virus}")
    else:
        return int(m.group('year'))

annotated_neut_titers = (
    neut_titers
    .assign(virus_year=lambda x: x['virus'].map(get_virus_year))
    )

annotated_neut_titers = (
    annotated_neut_titers
    .merge(annotated_neut_titers.query('virus_year == 1984')
                                [['serum', 'neut_titer']]
                                .drop_duplicates()
                                .rename(columns={'neut_titer': 'titer_1984'}),
           how='left', on='serum'
           )
    )
```

Get lower bound for neutralization titers and make sure consistent across all samples:


```python
titer_lower_bound = neut_titers.query('is_upper_bound')['neut_titer'].min()

print(f"The lower bound on measurable neutralization titers is {titer_lower_bound}")

assert (neut_titers.query('is_upper_bound')['neut_titer'] == titer_lower_bound).all()
```

    The lower bound on measurable neutralization titers is 10.0


Now plot neutralization titers as a function of virus isolation date for all viruses with a titer against the 1984 virus of at least 100.


```python
df = annotated_neut_titers.query('titer_1984 > 100')

p = (
    ggplot(df) +
    aes(x='virus_year', y='neut_titer') +
    geom_point() +
    geom_line() +
    geom_vline(data=df,
               mapping=aes(xintercept='collection_year'),
               color='orange',
               linetype='dashed',
               ) +
    facet_wrap('~ serum') +
    scale_y_log10() +
    theme_classic() +
    theme(axis_text_x=element_text(angle=90)) +
    geom_hline(yintercept=titer_lower_bound, linetype='dotted')
    )

_ = p.draw()
```


    
![png](analyze_neut_data_files/analyze_neut_data_23_0.png)
    



```python

```
