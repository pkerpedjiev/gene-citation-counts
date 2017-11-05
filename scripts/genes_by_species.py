import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="""
    python genes_by_species gene2pubmed

    Count how many citations there are for each species.
""")

    #parser.add_argument('argument', nargs=1)
    #parser.add_argument('-o', '--options', default='yo',
    #					 help="Some option", type='str')
    #parser.add_argument('-u', '--useless', action='store_true', 
    #					 help='Another useless option')
    parser.add_argument('filepath')

    args = parser.parse_args()
