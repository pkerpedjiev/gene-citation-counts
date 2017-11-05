#!/usr/bin/python

import collections as col
import operator
import sys
from optparse import OptionParser

def main():
    usage = """
    python count_by_column.py data_file column_number
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    if args[0] == '-':
        f = sys.stdin
    else:
        f = open(args[0], 'r')
    
    counts = col.defaultdict(int)

    for line in f:
        parts = line.split()
        counts[parts[int(args[1])-1]] += 1

    
    x = counts
    sorted_x = sorted(x.items(), key=operator.itemgetter(1))
    print "\n".join([" ".join(map(str,i)) for i in sorted_x])


if __name__ == '__main__':
    main()

