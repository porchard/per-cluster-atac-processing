#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('input_file')
parser.add_argument('--atac-barcodes', default='/home/porchard/github/snATACseq-NextFlow/737K-arc-v1.txt', help='default: /home/porchard/github/snATACseq-NextFlow/737K-arc-v1.txt')
parser.add_argument('--rna-barcodes', default='/home/porchard/github/snRNAseq-NextFlow/737K-arc-v1.txt', help='default: /home/porchard/github/snRNAseq-NextFlow/737K-arc-v1.txt')
parser.add_argument('--from', dest='fr', required=True, help='RNA or ATAC')
parser.add_argument('--to', required=True, help='RNA or ATAC')
parser.add_argument('--field', required=True, type=int, help='Column representing barcodes (indexed from 1)')
args = parser.parse_args()

#CLUSTERS = '/lab/work/porchard/daily/2020-12-21/labs/clusters.txt'

ATAC_BARCODE_LIST = args.atac_barcodes
RNA_BARCODE_LIST = args.rna_barcodes

assert(args.fr in ['ATAC', 'RNA'])
assert(args.to in ['ATAC', 'RNA'])
assert(args.to != args.fr)

# map RNA --> ATAC barcodes
atac_barcodes = pd.read_csv(ATAC_BARCODE_LIST, header=None, names=['atac_barcode'])
rna_barcodes = pd.read_csv(RNA_BARCODE_LIST, header=None, names=['rna_barcode'])
barcodes = pd.concat([atac_barcodes, rna_barcodes], axis=1)
conversions = dict(zip(barcodes.rna_barcode, barcodes.atac_barcode)) if (args.fr == 'RNA' and args.to == 'ATAC') else dict(zip(barcodes.atac_barcode, barcodes.rna_barcode))

input = pd.read_csv(args.input_file, sep='\t', header=None)
assert(len(input.columns) >= args.field)

# sometimes barcodes might be embedded in a larger string
# check if all elements in the column matches the conversion list
# if not, try to find an embedded barcode?? -- NOT YET IMPLEMENTED
if all(input[args.field - 1].isin(set(conversions.keys()))):
    input[args.field - 1] = input[args.field - 1].map(conversions)
else:
    raise ValueError('Not all elements in the field to convert are in the barcode list')
    #for i in input[args.field - 1]:


input.to_csv(sys.stdout, sep='\t', index=False, header=None)