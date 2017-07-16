#!/usr/bin/python

from __future__ import print_function

import gzip
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="""
    
    python add_gene_id_to_refseq.py human_refgene.txt.gz human_gene2refseq
""")

    parser.add_argument('gene2refseq')
    parser.add_argument('refGene')
    #parser.add_argument('-o', '--options', default='yo',
    #					 help="Some option", type='str')
    #parser.add_argument('-u', '--useless', action='store_true', 
    #					 help='Another useless option')

    args = parser.parse_args()
    refseq2gene = dict()

    with open(args.gene2refseq, 'r') as f:
        for line in f:
            parts = line.split('\t')
            gene_id = int(parts[1])
            refseq_id = parts[3].split('.')[0]

            refseq2gene[refseq_id] = gene_id

    with gzip.open(args.refGene, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[1] in refseq2gene:
                print("{}\t{}".format(str(refseq2gene[parts[1]]), "\t".join(parts)))
            else:
                print("not found:", parts[1], file=sys.stderr)
    

if __name__ == '__main__':
    main()


