#!/usr/bin/env python

import argparse
import pysam
import logging

parser = argparse.ArgumentParser()
parser.add_argument('bam_in', help = 'Bam file to adjust tags in')
parser.add_argument('bam_out_prefix', help = 'Output bam prefix (PREFIXcluster.bam)')
parser.add_argument('barcode_to_cluster_mapping', help = 'Barcode, cluster (tsv; barcodes should be in "CB" tag).')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


logging.info('Reading barcode to cluster mappings')
barcode_to_cluster = dict()
with open(args.barcode_to_cluster_mapping, 'r') as f:
    for line in f:
        barcode, cluster = line.rstrip().split()
        barcode_to_cluster[barcode] = cluster

# output new files
logging.info('Writing new bam files')
with pysam.AlignmentFile(args.bam_in, 'rb') as old_bam:
    bam_handles = dict()
    count = 0
    for read in old_bam.fetch(until_eof = True):
        count += 1
        if count % 1000000 == 0:
            logging.info('Processed {} reads'.format(count))
        if not read.has_tag('CB'):
            continue
        barcode = read.get_tag('CB')
        if barcode not in barcode_to_cluster:
            continue
        cluster = barcode_to_cluster[barcode]
        bam_name = '{}{}.bam'.format(args.bam_out_prefix, cluster)
        if bam_name not in bam_handles:
            bam_handles[bam_name] = pysam.AlignmentFile(bam_name, 'wb', template=old_bam)
        bam_handles[bam_name].write(read)

for nm, bam in bam_handles.items():
    bam.close()
            
logging.info('Done')
