"""Interactive multiple-sequence alignment viewer.

Code taken from:
https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner

"""

import collections

import os, io, random
import string
import numpy as np

from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

import panel as pn
import panel.widgets as pnw
pn.extension()

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot


def get_entropy(seqs):
    """Entropy for each site in list of sequences."""
    seqs = list(seqs)
    entropy = []
    for i in range(len(seqs[0])):
        counts = collections.Counter(s[i] for s in seqs)
        freqs = [n / len(seqs) for n in counts.values()]
        if len(freqs) == 1:
            entropy.append(0)
        else:
            entropy.append(-sum(np.array(freqs) * np.log(freqs)))
    return entropy


def get_colors(seqs):
    """color based on site entropy."""

    small_color = '#f76ab4'
    nucleophilic_color = '#ff7f00'
    hydrophobic_color = '#12ab0d'
    aromatic_color = '#84380b'
    acidic_color = '#e41a1c'
    amide_color = '#972aa8'
    basic_color = '#3c58e5'

    mapping_d = {'G':small_color, 'A':small_color,
                 'S':nucleophilic_color, 'T':nucleophilic_color, 'C':nucleophilic_color,
                 'V':hydrophobic_color, 'L':hydrophobic_color, 'I':hydrophobic_color, 'M':hydrophobic_color, 'P':hydrophobic_color,
                 'F':aromatic_color, 'Y':aromatic_color, 'W':aromatic_color,
                 'D':acidic_color, 'E':acidic_color,
                 'H':basic_color, 'K':basic_color, 'R':basic_color,
                 'N':amide_color, 'Q':amide_color,
                 '-':'gray'}

    entropy = get_entropy(seqs)
    colors = []
    for s in list(seqs):
        for i, a in enumerate(s):
            if entropy[i] > 0:
                colors.append(mapping_d[a])
            else:
                colors.append('white')
    return colors


def view_alignment(aln, fontsize="9pt", plot_width=800):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)    
    N = len(seqs[0])
    S = len(seqs)    
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    if N>100:
        viewlen=round(100 * plot_width / 800)
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"

    #entire sequence view (no text, with zoom)
    p = figure(title=None, plot_width= plot_width, plot_height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False  

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
                x_range=view_range, y_range=ids, tools="xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p],[p1]], toolbar_location='below')
    return p
