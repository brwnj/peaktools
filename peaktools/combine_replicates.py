#! /usr/bin/env python

'''peaktools-combine-replicates: combine replicate peak calls to generate
merged peak calls.

Assumes BED files contain peak calls from the same strand.
'''

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '$Revision: 436 $'

# Copyright 2011, Jay R. Hesselberth

import sys
from collections import defaultdict, Counter
from seqtools.formats import read_bed
from genomedata._util import maybe_gzip_open

def combine_replicates(bedfilenames, min_rep, verbose):

    # counts the number of times a base is covered by a peak call
    region_counts = defaultdict(Counter)

    for bedfilename in bedfilenames:
        if verbose:
            print >>sys.stderr, ">> loading regions from %s" % \
                bedfilename

        with maybe_gzip_open(bedfilename) as bedfile:
            for datum in read_bed(bedfile):
                for pos in range(datum.chromStart, datum.chromEnd):
                    region_counts[datum.chrom][pos] += 1

    # go over the merged counts and report contiguous regions above
    # min_rep
    region_start = region_stop = None

    for chrom in sorted(region_counts):
        for pos in sorted(region_counts[chrom]):

            count = region_counts[chrom][pos]

            if count >= min_rep:
                if not region_start:
                    region_start = pos
                region_stop = pos     

            elif count < min_rep and region_start:

                fields = (chrom, region_start, region_stop)
                print '\t'.join(map(str, fields))

                region_start = region_stop = None

def parse_options(args):
    from optparse import OptionParser
    
    description = ("Merge peak calls from BED regions")
    usage = '%prog [options] BEDFILES...'
    version = '%%prog %s' % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-m", "--minimum-replicated",
        action="store", type='int', dest="min_rep",
        help="minimum times a base is in a peak [default: %default]")

    parser.add_option("-v", "--verbose",
        action="store_true", dest="verbose",
        help="maximum verbosity [default: %default]")

    parser.set_defaults(min_rep=2,
                        verbose=False)

    options, args = parser.parse_args()

    if len(args) < 2:
        parser.error('specify at least 2 BEDFILES')

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    bedfilenames = args

    kwargs = {"min_rep":options.min_rep,
              "verbose":options.verbose}

    return combine_replicates(bedfilenames, **kwargs)

if __name__ == '__main__':
    sys.exit(main())

