#!/usr/bin/env python

import sys
import copy
import os
import json
import gzip
import numpy
import argparse
import re
import pysam
import random
import logging

parser = argparse.ArgumentParser()
parser.add_argument('bam_in', help = 'Bam file to adjust tags in')
parser.add_argument('bam_out', help = 'Output bam')
parser.add_argument('barcode_list', help = 'List of barcodes to keep (barcodes should be in "CB" tag).')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


logging.info('Reading list of barcodes to keep')
barcodes_keep = set()
with open(args.barcode_list, 'r') as f:
    for line in f:
        barcodes_keep.add(line.rstrip())


# output new file 
logging.info('Writing new bam file ({})'.format(args.bam_out))
with pysam.AlignmentFile(args.bam_in, 'rb') as old_bam:
    with pysam.AlignmentFile(args.bam_out, 'wb', template=old_bam) as new_bam:
        keep_count = 0
        remove_count = 0
        for read in old_bam.fetch(until_eof = True):
            if not read.has_tag('CB'):
                continue
            barcode = read.get_tag('CB')
            if (keep_count + remove_count) % 1000000 == 0:
                logging.info('Processed {} reads...'.format(keep_count + remove_count))
            if barcode in barcodes_keep:
                keep_count += 1
                new_bam.write(read)
            else:
                remove_count += 1

logging.info('Finished filtering {}. Kept {} reads and removed {}'.format(args.bam_in, keep_count, remove_count))
