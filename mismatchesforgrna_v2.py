## author Tanaya Jadhav
##compiles mismatches for grna sequences into a single tsv file
##needs manual changes for more than 2 mismatches
import itertools, random
from cacheout import Cache
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import glob
import os

all_c=set('AGCT')
get_other = lambda x : list(all_c.difference(x))
# print(all_c)
# print(get_other)

other={}
for c in all_c:
    other[c]=get_other(c)
# print(other)

def get_changed(sub, i):
    return [sub[0:i]+c+sub[i+1:] for c in other[sub[i]]]

mmatchHash=Cache(maxsize=256*256, ttl=0, timer=time.time, default=None)


def get_mismatch(d, setOfMmatch):
    if d == 0:
        return setOfMmatch

    newMmatches = []
    count = 0
    for sub in setOfMmatch:
        count+=1
        # print('setOfMmatch:', setOfMmatch, 'sub:',sub)
        newMmatches.extend(list(map(lambda x: ''.join(x),
                                    itertools.chain.from_iterable(([get_changed(sub, i) for i, c in enumerate(sub)])))))
        # print(newMmatches)
    if count == 1:
        mismatch1set = newMmatches
    else:
        mismatch2set = newMmatches
    setOfMmatch = setOfMmatch.union(newMmatches)
    if not mmatchHash.get((d - 1, str(setOfMmatch)), 0):
        mmatchHash.set((d - 1, str(setOfMmatch)), get_mismatch(d - 1, setOfMmatch))

    return (mmatchHash.get((d - 1, str(setOfMmatch))))


## load seqkit results
headercols = ['seqid', 'queryid', 'queryseq', 'strand', 'start', 'end', 'matched']
seqkitdf = pd.read_csv('cat_haroldReference_probes_0-500kb.txt', sep='\t', header=None, names=headercols)
# print(seqkitdf.head(5))
seqkitdf['twentymer'] = seqkitdf['matched'].str[0:20]
seqkitdf['twentymer'] = seqkitdf['twentymer'].str.upper()

## load probes
probedf = pd.read_csv('linked_probes.txt', sep='\t')
# print(probedf.head(5))
grnas = list(probedf['probe_without_NGG'])
# print(grnas)

## create output directory
outdir = 'linkedprobes_mismatchoutput/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

## loop through all grnas to generate mismatches
for grna in grnas:
    dna = grna.upper()
    for i in range(0,3):
        hamm_dist = i
        length = len(dna)
        mismatchlist = list(itertools.chain.from_iterable(
            [get_mismatch(hamm_dist, {dna[i:i + length]}) for i in range(len(dna) - length + 1)]))

        ##check against all 20merNGG from previous seqkit results
        subdf = seqkitdf[seqkitdf['twentymer'].isin(mismatchlist)]
        subdf.to_csv(outdir+grna+'_'+str(hamm_dist)+'mismatches.tsv', sep='\t', index=None)

    perfectmatchesdf = pd.read_csv(outdir+grna+'_0mismatches.tsv', sep='\t', header=0)
    perfectmatches = set(perfectmatchesdf['start'])
    # print(perfectmatches)
    mis1df = pd.read_csv(outdir+grna+'_1mismatches.tsv', sep='\t', header=0)
    mis1 = set(mis1df['start'])
    # print(mis1)
    truemis1 = list(mis1 - perfectmatches)
    mis2df = pd.read_csv(outdir+grna+'_2mismatches.tsv', sep='\t', header=0)
    mis2 = set(mis2df['start'])
    truemis2 = list(mis2 - mis1)
    mis2df['mismatches'] = ''
    mis2df.loc[mis2df.start.isin(perfectmatches), 'mismatches'] = 0
    mis2df.loc[mis2df.start.isin(truemis1), 'mismatches'] = 1
    mis2df.loc[mis2df.start.isin(truemis2), 'mismatches'] = 2
    mis2df.to_csv(outdir+grna+'_upto2mismatches.tsv', sep='\t', index=None)
    ##delete intermediate files
    if os.path.exists(outdir+grna+'_0mismatches.tsv'):
        os.remove(outdir+grna+'_0mismatches.tsv')
    if os.path.exists(outdir+grna+'_1mismatches.tsv'):
        os.remove(outdir+grna+'_1mismatches.tsv')
    if os.path.exists(outdir+grna+'_2mismatches.tsv'):
        os.remove(outdir+grna+'_2mismatches.tsv')

