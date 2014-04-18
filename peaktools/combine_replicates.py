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
import itertools
from toolshed import reader
from collections import defaultdict, Counter


def print_peak(chrom, start, stop):
    print "{chrom}\t{start}\t{stop}".format(chrom=chrom,
                                            start=start,
                                            stop=stop)


def get_region_counts(bedfilenames, verbose):
    # counts the number of times a base is covered by a peak call
    region_counts = defaultdict(Counter)

    for bedfilename in bedfilenames:
        if verbose:
            print >>sys.stderr, ">> loading regions from %s" % \
                bedfilename

        for d in reader(bedfilename, header=['chrom','start','end','name','score','strand']):
            for pos in range(int(d['start']), int(d['end'])):
                region_counts[d['chrom']][pos] += 1

    return region_counts


def group_ranges(count_items):
    # continuous intervals
    return (list(g) for k, g in itertools.groupby(count_items,\
                lambda x, y=itertools.count(): x[0]-next(y)))


def combine_replicates(bedfilenames, min_rep, verbose):

    region_counts = get_region_counts(bedfilenames, verbose)

    # go over the merged counts and report contiguous regions above
    # min_rep
    for chrom in sorted(region_counts):

        # grouped by continuous range
        for grouped in group_ranges(sorted(region_counts[chrom].items(),\
                                    key=lambda x: x[0])):

            region_start, region_stop, previous_pos = None, None, None
            peaks = []

            for pos, count in grouped:

                if count >= min_rep:
                    if not region_start:
                        region_start = pos
                    region_stop = pos + 1

                elif region_start:
                    print_peak(chrom, region_start, region_stop)
                    region_start, region_stop = None, None

            if region_start:
                print_peak(chrom, region_start, region_stop)


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
