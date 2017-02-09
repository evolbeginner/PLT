#! /bin/env python


import sys
import re

from Bio import SeqIO


#######################################################
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    if re.search(r'[^ATGCNatgcn]', str(seq_record.seq)):
        print seq_record.id



