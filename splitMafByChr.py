#! /bin/env python


import sys
import getopt
import re
import os

import selfbasics


#####################################################################
infile = None
outdir = None
is_force = False


#####################################################################
def splitMaf(infile, outdir):
    out_fhs = {}
    is_new_block = True
    counter = 0
    first_line = None
    in_fh = open(infile, 'r')
    for line in in_fh:
        if line.startswith('#'):
            continue
        line = line.rstrip('\n\r')
        counter += 1
        if counter == 1:
            first_line = line
        if counter == 2:
            #s ATH.2           558921 442 + 19698289 ATGGAGCCTTTA
            fields = re.split("\s+", line)
            chr_full_name = fields[1]
            m = re.search("([^.]+)\.([^.]+)", chr_full_name) 
            chr = m.group(2)
            outfile = os.path.join(outdir, chr + "." + 'maf')
            out_fh = open(outfile, 'a')
            out_fhs[chr] = out_fh
            out_fh.write(first_line+'\n')
            out_fh.write(line+'\n')
        elif counter > 2:
            out_fh.write(line+'\n')
            if line == '':
                counter = 0
    for chr,h in out_fhs.items():
        h.close()
    in_fh.close()


#####################################################################
try:
    opts, args = getopt.getopt(
        sys.argv[1:],
        'i:',
        ['in=','outdir=','force'],
    )
except getopt.GetoptError:
    print "Illegal arguments!"
    print "Exiting ......"


for opt, value in opts:
    if opt == '-i' or opt == '--in':
        infile = value
    elif opt == '--outdir':
        outdir = value
    elif opt == '--force':
        is_force = True
    else:
        print opt
        print "Illlegal argument! Exiting ......"
        sys.exit()


if not outdir:
    print "outdir has to be given! Exiting ......"
    sys.exit()


selfbasics.make_dir_force(outdir, is_force)


#####################################################################
splitMaf(infile, outdir)


