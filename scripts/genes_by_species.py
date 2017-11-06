import argparse
import collections as col
import sys

def main():
    parser = argparse.ArgumentParser(description="""
    python genes_by_species gene2pubmed

    Count how many citations there are for each species.
""")

    parser.add_argument('filename')
    #parser.add_argument('-o', '--options', default='yo',
    #					 help="Some option", type='str')
    #parser.add_argument('-u', '--useless', action='store_true', 
    #					 help='Another useless option')

    args = parser.parse_args()

    if args.filename == '-':
        f = sys.stdin
    else:
        f = open(args.filename)

    citations_by_taxid = col.defaultdict(set)

    for line in f:
        # usually the top line will contain header information
        if line[0] == '#':
            continue
        parts = line.strip().split()
        taxid = parts[0]
        citation = parts[2]

        citations_by_taxid[taxid].add(citation)

    for key in citations_by_taxid:
        print("{}\t{}".format(key, len(citations_by_taxid[key])))

if __name__ == "__main__":
    main()
