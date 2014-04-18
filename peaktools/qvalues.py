#! /usr/bin/env python

'''peaktools-qvalues: Calculate q-values from BED files
containing peak and null p-values.

Assumes that p-values are stored in the BED score field.

See Chen et al. (2010) Bioinformatics 26:i334-i342 for a description of
the algorithm.'''

import sys

from numpy import isnan, log10
from toolshed import reader
from segtools import ProgressBar

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = 'Revision: XXX'

# Copyright 2010,2011 Jay R. Hesselberth


def calc_qvalues(real_bedfilename, null_bedfilename, log_pvalues, verbose):

    # read in real p-values.
    real_pvals = read_pvalues(real_bedfilename, log_pvalues, verbose)

    # read in null p-values
    null_pvals = read_pvalues(null_bedfilename, log_pvalues, verbose)

    num_real = float(len(real_pvals))
    num_null = float(len(null_pvals))

    # make sure both are defined
    assert num_real and num_null

    # normalization factor to account for different numbers of real and
    # null p-values
    frac_real = num_real / num_null

    if verbose:
        print >>sys.stderr, ">> normalization factor: %.5f" % frac_real

    # compute pvalue thresholds
    pval_thresh = compute_pval_thresh(real_pvals, null_pvals, verbose)

    # go back over real pvalues and assign qvalues
    for d in reader(real_bedfilename, header=['chrom','start','end','name','score','strand']):

        if log_pvalues:
            pval = float(d['score'])
        else:
            pval = -1 * log10(pval)

        qval = pval_thresh[pval]

        norm_qval = qval * frac_real

        # print in table format
        fields = (d['chrom'], d['start'], d['end'],
                  d['name'], pval, d['strand'], norm_qval)
        print '\t'.join([str(f) for f in fields])


def read_pvalues(bedfilename, log_pvalues, verbose):
    ''' read in p-values from a bed file score field.

    returns: list sorted by signifance (most significant first)'''
    pvals = []

    if verbose:
        print >>sys.stderr, ">> reading p-values from %s .." % bedfilename

    for d in reader(bedfilename, header=['chrom','start','end','name','score','strand']):
        if log_pvalues:
            pval = float(d['score'])
        else:
            pval = -1 * log10(pval)
        pvals.append(pval)

    if verbose:
        print >>sys.stderr, ">> read %d p-values" % len(pvals)

    # sort the pvalues from most to least signif (smallest to largest) and
    # reverse so largest are first
    pvals.sort()

    # if pvals are log transformed, biggest (i.e. most significant) are
    # first
    if log_pvalues: pvals.reverse()

    return pvals


def compute_pval_thresh(real_pvals, null_pvals, verbose):
    ''' given an observed p-value and a collection of real and random p-values,
    calculate the empirical false discovery rate (FDR) for the p-value.

    returns: dictionary of floats'''

    pval_thresh = dict()

    if verbose:
        label = ">> computing thresholds: "
        progress = ProgressBar(len(real_pvals), label=label)

    for obs_pval in real_pvals:

        if verbose: progress.next()

        if obs_pval in pval_thresh: continue

        real_pvals_better = num_pvals_better(obs_pval, real_pvals, verbose)
        null_pvals_better = num_pvals_better(obs_pval, null_pvals, verbose)

        if real_pvals_better == 0:
            qval = 0.0
        else:
            qval = null_pvals_better / real_pvals_better

        pval_thresh[obs_pval] = qval

    if verbose: progress.end()

    return pval_thresh


def num_pvals_better(obs_pval, pvals, verbose):
    ''' calculate the number of p-values in a collection that are
    more significant than an observed p-value.

    obs_pval: observed pvalue
    pvals: sorted list of pvalues
    verbose: verbosity flag

    returns: float'''

    num_better = 0.0

    for pval in pvals:
        if pval > obs_pval:
            num_better += 1
        else:
            break # can't get any better

    return num_better


def parse_options(args):
    from optparse import OptionParser

    description = ("Compute q-values (i.e. FDR's) from real and null p-values "
                   "in BED format")
    usage = '%prog [options] REAL_BEDFILENAME NULL_BEDFILENAME'
    version = '%%prog %s' % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    # XXX this option makes no sense
    parser.add_option("-l", "--log-pvalues",
        action="store_true", dest="log_pvalues",
        help="pvalues are log-transformed [default: %default]",
        default=True)

    parser.add_option("-v", "--verbose",
        action="store_true", dest="verbose",
        help="verbose output [default: %default]",
        default=False)

    options, args = parser.parse_args()

    if len(args) < 2:
        parser.error('specify REAL_BEDFILENAME and NULL_BEDFILENAME')

    return options, args


def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    real_bedfilename, null_bedfilename = args
    kwargs = {'log_pvalues':options.log_pvalues,
              'verbose':options.verbose}
    return calc_qvalues(real_bedfilename, null_bedfilename, **kwargs)


if __name__ == '__main__':
    sys.exit(main())
