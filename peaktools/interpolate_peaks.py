#! /usr/bin/env python

'''peaktools-interpolate-peaks: XXX
'''

import pdb
import sys
import warnings

from numpy import nan_to_num
from segway.bed import read_native
from genomedata import Genome
from genomedata._util import maybe_gzip_open

from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

rget = robjects.r('get')
renv = robjects.r('environment')
stats = importr('stats')
spline = stats.spline
splinefun = stats.splinefun

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '$Revision: 439 $'

# Copyright 2010,2011 Jay R. Hesselberth

def interpolate_peaks(gdfilename, bedfilename, trackname, spec_chrom,
                      verbose):

    warnings.simplefilter('ignore')
    with Genome(gdfilename) as genome, \
        maybe_gzip_open(bedfilename) as bedfile:
        for datum in read_native(bedfile):
        
            chrom = datum.chrom
            peak_start = datum.chromStart
            peak_end = datum.chromEnd

            # can parallelize by chrom
            if spec_chrom and spec_chrom != chrom: continue

            xs = range(peak_start, peak_end)
            ys = genome[chrom][peak_start:peak_end, trackname]

            interp_start, interp_end = fit_spline(xs, ys)

            fields = (chrom, interp_start, interp_end,
                      datum.name, datum.score, datum.strand)
            print '\t'.join(map(str, fields))

def fit_spline(xs, ys):
    ''' fit a spline curve to peak signal. XXX'''

    # coerce the y values to list of ints
    ys = list(nan_to_num(ys).astype(int))

    result = splinefun(xs, ys)
    # adapted from the R man page on splinefun
    coefs = rget("z", envir=renv(result))

    coefs_b = [i for i in coefs[3]]
    coefs_c = [i for i in coefs[4]]
    coefs_d = [i for i in coefs[5]]

    b_max_idx = max([(coef, idx) for idx, coef in 
                  enumerate(coefs_b)])[1]
    c_max_idx = max([(coef, idx) for idx, coef in
                  enumerate(coefs_c)])[1]

    # print ">> c max: %d -> b max: %d" % (c_max, b_max)
    # pdb.set_trace()
   
    return (xs[c_max_idx], xs[b_max_idx])

def parse_options(args):
    from optparse import OptionParser, OptionGroup
    
    description = ("Find peak changepoints by interpolation.")
    usage = '%prog [options] GENOMEDATADIR BEDFILENAME'
    version = '%%prog %s' % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Required")
    group.add_option("-t", "--trackname",
        action="store", dest="trackname",
        help="trackname [default: %default]",
        default=None)    

    parser.add_option_group(group)

    group = OptionGroup(parser, "Variables")
    group.add_option("-c", "--chrom",
        action="store", type='string', dest="spec_chrom",
        help="analyze specific chromosome [default: %default]",
        default=None)

    parser.add_option_group(group)

    group = OptionGroup(parser, "Flags")
    group.add_option("-v", "--verbose",
        action="store_true", dest="verbose",
        help="maximum verbosity [default: %default]",
        default=False)    

    parser.add_option_group(group)

    options, args = parser.parse_args()

    if len(args) != 2:
        parser.error('specify GENOMEDATADIR and BEDFILENAME')
    if not options.trackname:
        parser.error('specify trackname')

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    gdfilename, bedfilename = args
    kwargs = {'trackname':options.trackname,
              'spec_chrom':options.spec_chrom,
              'verbose':options.verbose}

    return interpolate_peaks(gdfilename, bedfilename, **kwargs)

if __name__ == '__main__':
    sys.exit(main())

