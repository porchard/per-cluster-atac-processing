#!/usr/bin/env python
# coding: utf-8

import os
import sys
import matplotlib.pyplot as plt
import pygenometracks.tracks as pygtk
import argparse
import glob
import pandas as pd
import re
import pybedtools as bt

parser = argparse.ArgumentParser()
parser.add_argument('--gene-bed', default=None, dest='gene_bed', help='Path to GTF/gene bed file. Used for gene models.')
parser.add_argument('--out', default='marker-genes-atac.png', help='Name of output plot.')
parser.add_argument('--cluster-names', dest='cluster_names', default=None)
parser.add_argument('--cluster-colors', dest='cluster_colors', default=None)
parser.add_argument('--exclude-clusters', dest='exclude_clusters', help='Cluster IDs to exclude', nargs='*', default=[])
parser.add_argument('--genes', nargs='+', help='List of genes to plot (in order)')
parser.add_argument('--bigwigs', nargs='+', help='Paths to bigwig files.')
args = parser.parse_args()

GENES = args.genes
BIGWIGS = args.bigwigs
GENE_BED = args.gene_bed
EXCLUDE_CLUSTERS = args.exclude_clusters


def plot_signal_at_genes(tracks, genes, gene_locations, gene_bed, height=None, width=None):
    """
    tracks should be a list of dicts. One element per bigwig, keys: 'name', 'color' (optional), 'bigwig'
    genes should be a list of genes to plot.
    gene_locations: should be a list of gene locations (chrom:start-end). If None, then infer from gene_bed
    gene_bed: gene bed file for drawing gene models.

    E.g.:
    #tracks = [{'name': 'Type II fibers', 'color': 'orange', 'bigwig': '/lab/work/porchard/sn-muscle-project/work/process-by-cluster/results/bigwigs/atac/0-hg19.bw'},
    #          {'name': 'Type I fibers', 'color': 'blue', 'bigwig': '/lab/work/porchard/sn-muscle-project/work/process-by-cluster/results/bigwigs/atac/1-hg19.bw'}]
    #genes = ['GAPDH', 'MYH1', 'MYH7']
    #gene_bed = '/lab/work/porchard/sn-muscle-project/data/gencode-bed/gencode.hg19.genes.bed'
    #gene_locations = {'GAPDH': 'chr12:6642596-6648525', 'MYH1': 'chr17:10414886-10423629', 'MYH7': 'chr14:23878228-23912613'}
    """

    tk = None
    if gene_bed is not None:
        track_config = dict(file=gene_bed, file_type='gtf', color='black', merge_transcripts=True, labels=False, display='interleaved',arrow_interval=999999, section_name='', prefered_name='gene_name', font_size=8, style="UCSC", all_labels_inside=True, labels_in_margin=False)
        tk = pygtk.BedTrack(track_config)

    # lay out figure
    # columns: first column is the label, then one column per gene
    # rows: first row is the gene name, second row in the gene model, then one row per cell type
    HEIGHT_RATIOS = [0.5, 1] + [3 for i in range(len(tracks))]
    fig, axs = plt.subplots(ncols=1+len(genes), nrows=2+len(tracks), gridspec_kw={'hspace':0, 'wspace':0, 'height_ratios': HEIGHT_RATIOS}, figsize=(len(genes)*1.2 if width is None else width, len(tracks)*0.6 if height is None else height))
    axs[0,0].axis('off')
    axs[1,0].axis('off')
    for ax in axs[:,1]:
        ax.spines['left'].set_visible(False)
    for ax in axs[1,:]:
        ax.spines['bottom'].set_visible(False)
    for ax in axs[2,:]:
        ax.spines['top'].set_visible(False)

    # fill in the gene names
    for gene, ax in zip(genes, axs[0:,1:].flatten()):
        ax.text(x=0.5, y=0, s=gene, ha='center', va='bottom', fontsize='large', fontstyle='italic')
        for t in ax.get_yticklabels():
            t.set(visible=False)
        ax.axis('off')


    # fill in the gene models
    for gene, ax in zip(genes, axs[1:,1:].flatten()):
        chrom, start, end = re.match('^(.*):(\d+)-(\d+)$', gene_locations[gene]).groups()
        if tk is not None:
            tk.plot(ax, chrom, int(start), int(end))
        for t in ax.get_yticklabels():
            t.set(visible=False)
        ax.axis('off')


    # fill in cell type names
    for t, ax in zip(tracks, axs[2:,0].flatten()):
        ax.text(x=1, y=0.5, s=t['name'], ha='right', va='center', fontsize='large')
        ax.axis('off')
    

    # fill in the bigwig tracks
    for track_index, t in enumerate(tracks, 2):
        for gene_index, gene in enumerate(genes, 1):
            properties_dict = {'file': t['bigwig'], 'height': 3, 'color': t['color'] if 'color' in t else 'black'}
            bw = pygtk.BigWigTrack(properties_dict)
            chrom, start, end = re.match('^(.*):(\d+)-(\d+)$', gene_locations[gene]).groups()
            ax = axs[track_index,gene_index]
            bw.plot(ax, chrom, int(start), int(end))
            ax.xaxis.set_ticklabels([])
            ax.set_yticks([], [])
            ax.set_xticks([], [])
            ax.set_ylim(0, 6)
            ax.yaxis.set_ticklabels([])

    fig.tight_layout()
    return (fig, axs)




tracks = [{'name': 'cluster_' + os.path.basename(i).replace('.bw', ''), 'bigwig': i} for i in BIGWIGS]
tracks = [i for i in tracks if i['name'] not in [f'cluster_{j}' for j in EXCLUDE_CLUSTERS]]
if args.cluster_names is not None:
    cluster_names = pd.read_csv(args.cluster_names, sep='\t', dtype=str, names=['cluster_id', 'cluster_name']).set_index('cluster_id').cluster_name.to_dict()
    cluster_colors = pd.read_csv(args.cluster_colors, sep='\t', dtype=str, names=['cluster_id', 'cluster_color']).set_index('cluster_id').cluster_color.to_dict()
    for i in range(len(tracks)):
        cluster_id = tracks[i]['name'].replace('cluster_', '')
        tracks[i]['name'] = '{}'.format(cluster_names[cluster_id])
        tracks[i]['color'] = '{}'.format(cluster_colors[cluster_id])


tmp = pd.read_csv(GENE_BED, sep='\t', header=None).iloc[:,[0, 1, 2, 3]]
tmp.columns = ['chrom', 'start', 'end', 'gene']

gene_locations = dict()
for gene, df in tmp[tmp.gene.isin(GENES)].groupby('gene'):
    chrom, start, end = bt.BedTool().from_dataframe(df).sort().merge().to_dataframe().iloc[0,:].to_list()
    start = int(start) - 5000
    end = int(end) + 5000
    gene_locations[gene] = f'{chrom}:{start}-{end}'


fig, axs = plot_signal_at_genes(tracks, GENES, gene_locations, GENE_BED)
fig.savefig(args.out, facecolor='white')
fig.clf()
