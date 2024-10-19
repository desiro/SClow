#!/usr/bin/env python3
# script: levenshtein.py
# author: Daniel Desiro'

import Levenshtein as ls
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import sys
import os
import re

infasta = sys.argv[1]
outfasta = sys.argv[2]
addfasta = sys.argv[3]


def main():
    addreads = readFasta(addfasta)
    inreads = readFasta(infasta)
    reads = {k.split("|")[1]:v for i,(k,v) in enumerate(inreads.items())}
    reads = {k:v for i,(k,v) in enumerate(reads.items()) if i < 3} | addreads | {k:v for i,(k,v) in enumerate(reads.items()) if i >= 3}
    heatMat = np.zeros(shape=(len(reads),len(reads)))
    for i,(k1,s1) in enumerate(reads.items()):
        for j,(k2,s2) in enumerate(reads.items()):
            heatMat[i,j] = ls.distance(s1,s2)
    s = 1.0
    n, m = heatMat.shape
    mx = (4.8 / n) * m + 1.6
    fig, ax = plt.subplots(figsize=(mx*s,4.8*s))
    im, cbar = heatmap(heatMat, reads.keys(), reads.keys(), ax=ax, cmap="YlOrRd", cbarlabel="Levenshtein distance")
    texts = annotate_heatmap(im, valfmt="{x:.0f}")
    fig.tight_layout()
    fig.savefig(f"{outfasta}-Levenshtein_heatmap.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    fig.savefig(f"{outfasta}-Levenshtein_heatmap.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(fig)




def readFasta(infasta):
    fasta_dict, RNA = dict(), ""
    with open(infasta, "r") as infa:
        for line in infa:
            line = line.strip()
            if re.match(r">", line):
                if RNA: fasta_dict[name] = revComp(RNA)
                RNA, name = "", line[1:].split()[0]
            else:
                RNA += line
        fasta_dict[name] = revComp(RNA)
    return fasta_dict

def revComp(RNA, cmp=False, rev=False):
    RNA = RNA.upper()
    D2Rc = {"A":"U","T":"A","U":"A","C":"G","G":"C","R":"Y","Y":"R","M":"K",\
            "K":"M","S":"W","W":"S","B":"V","V":"B","D":"H","H":"D","N":"N"}
    if cmp: RNA = "".join(D2Rc[i] for i in RNA)
    else:   RNA = RNA.replace("T","U")
    if rev: RNA = RNA[::-1]
    return RNA

def heatmap(data, row_labels, col_labels, ax=None, cbar_kw=None, cbarlabel="", **kwargs):
    if ax is None:      ax = plt.gca()
    if cbar_kw is None: cbar_kw = {}
    im = ax.imshow(data, **kwargs)
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")
    ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}", textcolors=("black", "white"), threshold=None, **textkw):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()
    if threshold is not None: threshold = im.norm(threshold)
    else:                     threshold = im.norm(data.max())/2.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)
    if isinstance(valfmt, str):
        valfmt = StrMethodFormatter(valfmt)
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)
    return texts

main()