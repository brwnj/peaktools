#! /usr/bin/env python

'''peaktools-qvalue-summary: print summary report for qvalue
calculations'''

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '$Revision: 445 $'

# Copyright 2010,2011 Jay R. Hesselberth

import sys
from tabdelim import DictReader
from genomedata._util import maybe_gzip_open

DEFAULT_THRESH = [0.0, 0.001, 0.01, 0.05, 0.1, 0.25,
                 0.5, 0.75, 0.9, 1.0]

QVAL_FIELDNAMES = ['chrom', 'chromStart', 'chromEnd',
                   'trackname', 'pvalue', 'strand', 'qvalue']

def qvalue_summary(qvalue_filename, thresholds):

    qvalues = read_qvalues(qvalue_filename)

    print '# total features: %d' % len(qvalues)
    print '\t'.join(['# qvalue_threshold','num_features'])

    for threshold in sorted(thresholds):
        num_better = num_better_qvalues(qvalues, threshold)
        print '\t'.join(map(str,[threshold,num_better]))

def num_better_qvalues(qvalues, threshold):
    ''' calculate the number of qvalues better than or equal to a threshold.

    returns: int'''
    num_better = 0
    for qval in qvalues:
        if qval <= threshold:
            num_better += 1
        else:
            break # qvalues are sorted so can do this safely
    return num_better

def read_qvalues(qvalue_filename):
    ''' read in q-values from file.

    returns: sorted list of q-values'''
    qvalues = []

    with maybe_gzip_open(qvalue_filename) as qvalue_file:
        for datum in DictReader(qvalue_file, QVAL_FIELDNAMES):
            qvalues.append(float(datum['qvalue']))

    qvalues.sort()

    return qvalues

def parse_options(args):
    from optparse import OptionParser

    description = ("report a summary of calculate q-values")
    usage = '%prog [options] QVALUE_DATA_FILENAME'
    version = '%%prog %s' % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    thresholds = ', '.join(map(str, DEFAULT_THRESH))
    parser.add_option("-t", "--threshold", action="append", default=[],
        metavar="THRESHOLD",
        help="append THRESHOLD to list of thresholds to analyze"
        " (default: %s)" % thresholds)

    options, args = parser.parse_args()

    if len(args) < 1:
        parser.error('specify QVALUE file')

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    qvalue_filename = args[0]

    if len(options.threshold) == 0:
        thresholds = DEFAULT_THRESH
    else:
        thresholds = [float(i) for i in options.threshold]

    return qvalue_summary(qvalue_filename, thresholds)

if __name__ == '__main__':
    sys.exit(main())

