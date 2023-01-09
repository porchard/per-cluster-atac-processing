#!/usr/bin/env python
# coding: utf-8

import sys
import os
import pybedtools as bt
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Convert a set of summits to a set of peaks as described in methods of paper "Epigenomic State Transitions Characterize Tumor Progression in Mouse Lung Adenocarcinoma"', add_help=True)
parser.add_argument('summits_bed', type=str, help='Bed file of summits. Chrom, start, end, name (must be unique), -log10(q). As output from MACS2 callpeak --call-summits')
parser.add_argument('--extension', type=int, default=150, help='Extend summits by this amount on either side (default: 150)')
args = parser.parse_args()

SUMMITS = args.summits_bed
EXTENSION = args.extension

if not EXTENSION >= 0:
    raise ValueError('--extension must be a value greater than 0')

summits = pd.read_csv(SUMMITS, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score'])

# extend summits --> peaks
peaks = summits.copy()
peaks.start = peaks.start - EXTENSION
if peaks.start.min() < 0:
    sys.stderr.write('Truncating {} peaks that would have start sites < 0\n'.format(len(peaks[peaks.start<0])))
    peaks.start = peaks.start.map(lambda x: x if x >= 0 else 0)
peaks.end = peaks.end + EXTENSION
peaks = bt.BedTool().from_dataframe(peaks).sort().saveas()

# find overlapping peaks
overlaps = peaks.intersect(peaks, wa=True, wb=True).to_dataframe()
overlaps.columns = ['chrom_1', 'start_1', 'end_1', 'name_1', 'score_1', 'chrom_2', 'start_2', 'end_2', 'name_2', 'score_2']
overlaps = overlaps[overlaps.name_1 != overlaps.name_2]
overlaps.head()
overlapping = dict() # peak --> [overlapping_peak_1, overlapping_peak_2, ...]
for p1, p2 in zip(overlaps.name_1, overlaps.name_2):
    if p1 not in overlapping:
        overlapping[p1] = []
    overlapping[p1].append(p2)

# drop less-significant overlapping peaks
summits = summits.sort_values('score', ascending=False)
drop_peaks = set() # drop these peaks in the end, as they overlap with peaks from more significant summits

for peak in summits.name:
    if peak in drop_peaks: # this peak overlaps a more significant peak -- so ignore it
        continue
    if peak not in overlapping: # there are no peaks overlapping this one, so keep it
        continue
    else:
        # each of the overlapping peaks need to be marked as DROP
        for overlapping_peak in overlapping[peak]:
            drop_peaks.add(overlapping_peak)

sys.stderr.write('Dropping {} of {} peaks\n'.format(len(drop_peaks), len(summits)))
summits.start = summits.start - EXTENSION
summits.end = summits.end + EXTENSION
if summits.start.min() < 0:
    sys.stderr.write('Truncating {} peaks that would have start sites < 0\n'.format(len(summits[summits.start<0])))
    summits.start = summits.start.map(lambda x: x if x >= 0 else 0)
summits[~summits.name.isin(drop_peaks)].to_csv(sys.stdout, sep='\t', index=False, header=None)

