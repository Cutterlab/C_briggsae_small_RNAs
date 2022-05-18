#!/usr/bin/env python2.7

# This script was written by Santiago Sánchez-Ramírez

# genSpecRef.py

import argparse
import re
import sys
import os
try:
    import pysam
except:
    sys.exit("pysam is not installed. Try: pip install pysam")

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="genSpecRef.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
Generates a genotype specific reference provided 
a VCF file and FASTA reference.""",
        epilog="""
All files must be indexed. So before running the code make sure 
that your reference FASTA file is indexed:

samtools faidx genome.fas

BGZIP compress and TABIX index your VCF file:

bgzip variants.vcf
tabix variants.vcf.gz

examples:
python genSpecRef.py -f genome.fas -v variants.vcf.gz -o my_output.fa

""")
    parser.add_argument(
    '--fasta', '-f', metavar='GENOME', type=str,
    help='FASTA file with the reference genome.')
    parser.add_argument(
    '--vcf', '-v', metavar='VCF', type=str, required=True,
    help='a tabix-indexed VCF file.')
    parser.add_argument(
    '--output', '-o', metavar='FILE', type=str,
    help='output file name.')
    #parser.add_argument(
    #'--blend', '-b', action="store_true", default=False,
    #help='concatenate GFF entries of FEAT into a single alignment. Useful for CDS. (default: False)')
    args = parser.parse_args()

    # read variant file and get samples
    vcf = pysam.VariantFile(args.vcf)
    # read genome reference file
    ref = pysam.FastaFile(args.fasta)
    # get fasta headers
    head = ref.references
    # open filehandle for output
    o = open(args.output, "w")
    # loop through reference headers
    for chrom in head:
        s = ref.fetch(chrom).upper()
        addposcum = 0
        for rec in vcf.fetch(chrom):
            alleles,addpos,refpos = UpdateAllele(rec)
            s = UpdateSeq(rec, alleles[1], addposcum, refpos, s)
            addposcum += addpos
            sys.stdout.write(" {}: {} \r".format(chrom,rec.pos)),
            sys.stdout.flush()
        o.write(">"+chrom+"1\n"+s+"\n")
        print ""
    o.close()

def UpdateAllele(vcfrec):
    allelelen = [ len(x) for x in vcfrec.alleles ]
    maxlen = allelelen.index(max(allelelen))
    alleles2 = tuple([ x + '-' * (allelelen[maxlen]-len(x)) for x in vcfrec.alleles ])
    return alleles2, abs(allelelen[0]-allelelen[maxlen]), allelelen[0]-1

def UpdateSeq(vcfrec, allele, addposcum, refpos, seq):
    pos = vcfrec.pos-1+addposcum
    seq = seq[:pos]+allele+seq[pos+1+refpos:]
    return seq

if __name__ == "__main__":
    main()

