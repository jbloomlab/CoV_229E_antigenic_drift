# Plot viral titers
This Python Jupyter notebook analyzes the neutralization data.

Import Python modules:


```python
import xml.etree.ElementTree as ElementTree

from IPython.display import display, HTML, SVG

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import svgutils.compose
```

Read in the titer data:


```python
titers = pd.read_csv('viral_titers.csv')

display(HTML(titers.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>comparison</th>
      <th>titer</th>
      <th>replicate</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1984</td>
      <td>virus</td>
      <td>4.401059e+07</td>
      <td>1</td>
    </tr>
    <tr>
      <td>1992</td>
      <td>virus</td>
      <td>3.674175e+07</td>
      <td>1</td>
    </tr>
    <tr>
      <td>2001</td>
      <td>virus</td>
      <td>3.252061e+07</td>
      <td>1</td>
    </tr>
    <tr>
      <td>2008</td>
      <td>virus</td>
      <td>1.521004e+07</td>
      <td>1</td>
    </tr>
    <tr>
      <td>2016</td>
      <td>virus</td>
      <td>3.721168e+07</td>
      <td>1</td>
    </tr>
    <tr>
      <td>RBD-chimera-1992</td>
      <td>virus</td>
      <td>8.596787e+06</td>
      <td>1</td>
    </tr>
    <tr>
      <td>RBD-chimera-2001</td>
      <td>virus</td>
      <td>3.018147e+07</td>
      <td>1</td>
    </tr>
    <tr>
      <td>RBD-chimera-2008</td>
      <td>virus</td>
      <td>2.807571e+07</td>
      <td>1</td>
    </tr>
    <tr>
      <td>RBD-chimera-2016</td>
      <td>virus</td>
      <td>2.568473e+07</td>
      <td>1</td>
    </tr>
    <tr>
      <td>APN only</td>
      <td>cells</td>
      <td>1.050000e-01</td>
      <td>1</td>
    </tr>
    <tr>
      <td>APN only</td>
      <td>cells</td>
      <td>9.600000e-02</td>
      <td>2</td>
    </tr>
    <tr>
      <td>TMPRSS2 only</td>
      <td>cells</td>
      <td>2.000000e-04</td>
      <td>1</td>
    </tr>
    <tr>
      <td>TMPRSS2 only</td>
      <td>cells</td>
      <td>1.000000e-05</td>
      <td>2</td>
    </tr>
    <tr>
      <td>APN + TMPRSS2</td>
      <td>cells</td>
      <td>1.000000e+00</td>
      <td>1</td>
    </tr>
    <tr>
      <td>APN + TMPRSS2</td>
      <td>cells</td>
      <td>1.000000e+00</td>
      <td>2</td>
    </tr>
    <tr>
      <td>full tail</td>
      <td>ZsGreen_titer</td>
      <td>4.000000e+03</td>
      <td>1</td>
    </tr>
    <tr>
      <td>full tail</td>
      <td>ZsGreen_titer</td>
      <td>1.700000e+03</td>
      <td>2</td>
    </tr>
    <tr>
      <td>tail deletion</td>
      <td>ZsGreen_titer</td>
      <td>1.100000e+04</td>
      <td>1</td>
    </tr>
    <tr>
      <td>tail deletion</td>
      <td>ZsGreen_titer</td>
      <td>1.000000e+04</td>
      <td>2</td>
    </tr>
    <tr>
      <td>no spike</td>
      <td>ZsGreen_titer</td>
      <td>0.000000e+00</td>
      <td>1</td>
    </tr>
    <tr>
      <td>no spike</td>
      <td>ZsGreen_titer</td>
      <td>0.000000e+00</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


Plot the titers:


```python
svgs = []

for group, xlabel, ylabel, ymin in [                         
                                    ('ZsGreen_titer', 'spike tail variant', 'transduction units / ml', 1e2),
                                    ('cells', 'cells expressing', 'relative titer', None),
                                    ('virus', 'spike', 'RLU / ml', 1e5),
                                    ]:
    df = titers.query('comparison == @group')
    if ymin is None:
        ymin = df['titer'].min()
        lod = None
    elif df['titer'].min() < ymin:
        lod = ymin
    else:
        lod = None
    p = (ggplot(df) +
         aes('name', 'titer') +
         geom_jitter(size=2.5, alpha=0.6, height=0, width=0.2 * (df['replicate'].nunique() - 1),
                     random_state=4) +
         xlab(xlabel) +
         scale_y_log10(name=ylabel,
                       limits=(0.67 * max(ymin, min(ymin, df['titer'].min())),
                               1.5 * df['titer'].max())) +
         theme_classic() +
         theme(axis_text_x=element_text(angle=90, size=10),
               figure_size=(0.5 * df['name'].nunique(), 2.2)
               )
         )
    if lod:
        p = p + geom_hline(yintercept=lod, linetype='dotted', color='gray', alpha=0.5, size=1)
    svg = f"results/titers_{ylabel.replace(' ', '_').replace('/', '_')}.svg"
    svgs.append(svg)
    print(f"Saving plot to {svg}")
    p.save(svg, verbose=False)
    fig = p.draw()
    display(fig)
    plt.close(fig)
```

    Saving plot to results/titers_transduction_units___ml.svg


    /fh/fast/bloom_j/software/miniconda3/envs/CoV_229E_antigenic_drift/lib/python3.8/site-packages/pandas/core/series.py:726: RuntimeWarning: divide by zero encountered in log10
    /fh/fast/bloom_j/software/miniconda3/envs/CoV_229E_antigenic_drift/lib/python3.8/site-packages/pandas/core/series.py:726: RuntimeWarning: divide by zero encountered in log10



    
![png](plot_titers_files/plot_titers_5_2.png)
    


    Saving plot to results/titers_relative_titer.svg



    
![png](plot_titers_files/plot_titers_5_4.png)
    


    Saving plot to results/titers_RLU___ml.svg



    
![png](plot_titers_files/plot_titers_5_6.png)
    


Use [svgutils](https://svgutils.readthedocs.io/) to assemble panels into a plot:


```python
titers_fig_svg = 'results/viral_titers_fig.svg'
print(f"Creating titer figure SVG as {titers_fig_svg}")

def svg_dim(svgfile, dim):
    """Get width or height `dim` of `svgfile` in points."""
    return float(ElementTree.parse(svgfile)
                            .getroot().attrib[dim]
                            .replace('px', '')
                            .replace('pt', '')
                            )

top_padding = 10
horiz_padding = 30

svgutils.compose.Figure(
    sum([svg_dim(svg, 'width') for svg in svgs]) + horiz_padding * (len(svgs) - 1),
    max([svg_dim(svg, 'height') for svg in svgs]) + top_padding,
    svgutils.compose.Panel(
        svgutils.compose.SVG(svgs[0]),
        svgutils.compose.Text('A', 2, 4, size=18, font='Arial'),
        ).move(0, top_padding),
    svgutils.compose.Panel(
        svgutils.compose.SVG(svgs[1]),
        svgutils.compose.Text('B', 2, 4, size=18, font='Arial'),
        ).move(svg_dim(svgs[0], 'width') + horiz_padding, top_padding),
    svgutils.compose.Panel(
        svgutils.compose.SVG(svgs[2]),
        svgutils.compose.Text('C', 2, 4, size=18, font='Arial'),
        ).move(svg_dim(svgs[0], 'width') + svg_dim(svgs[1], 'width') + 2 * horiz_padding,
               top_padding)
    ).save(titers_fig_svg)

display(SVG(titers_fig_svg))
```

    Creating titer figure SVG as results/viral_titers_fig.svg



    
![svg](plot_titers_files/plot_titers_7_1.svg)
    



```python

```
