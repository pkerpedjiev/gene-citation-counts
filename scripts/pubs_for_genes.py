#!/usr/bin/python

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="""
    
    python pubs_for_genes.py genes2pubmed genes_list

    List all the publications that reference one of the genes in 
    the genes_list argument
""")

    parser.add_argument('gene2pubmed')
    parser.add_argument('genes_list')

    #parser.add_argument('-o', '--options', default='yo',
    #					 help="Some option", type='str')
    #parser.add_argument('-u', '--useless', action='store_true', 
    #					 help='Another useless option')

    args = parser.parse_args()

    gene_list = set()

    with open(args.genes_list, 'r') as f:
        for line in f:
            gene_list.add(line.strip())

    #print("gene_list:", gene_list)

    all_human_papers = set()
    top100_gene_human_papers = set()

    with open(args.gene2pubmed, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if parts[0] == '9606':
                all_human_papers.add(parts[2])

            if parts[1] in gene_list:
                top100_gene_human_papers.add(parts[2])

    print("top100:", len(top100_gene_human_papers))
    print("all:", len(all_human_papers))

if __name__ == '__main__':
    main()


