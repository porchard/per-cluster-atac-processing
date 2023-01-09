#!/usr/bin/env python

import os
import sys
import re

# get transcript --> gene_name translations
GFF = sys.argv[1]
BED = sys.argv[2]

transcript_id_to_gene_name = dict()

with open(GFF, 'r') as f:
    for line in f:
        if re.match('^#', line):
            continue
        line = line.rstrip().split()
        attributes = [i.split('=') for i in line[-1].split(';')]
        for i in attributes:
            assert(len(i) == 2)
        attributes = {i[0]: i[1] for i in attributes}
        if 'transcript_id' not in attributes:
            continue
        transcript_id_to_gene_name[attributes['transcript_id']] = attributes['gene_name']

with open(BED, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        transcript_id = line[3]
        line[3] = transcript_id_to_gene_name[transcript_id] if transcript_id in transcript_id_to_gene_name else transcript_id
        print('\t'.join(line))
