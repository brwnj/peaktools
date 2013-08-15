#! /usr/bin/env python

'''peaktools-identify-peaks:  calculate poisson p-values for observations in a
genomedata track.  essentially implements previous poisson peak callers
e.g. MACS and Gene Yeo papers
'''

import pdb
import sys
import warnings
from itertools import combinations

from numpy import nansum, isnan
from numpy.random import shuffle
from genomedata import Genome
from segtools import ProgressBar

from rpy2.robjects.packages import importr
stats = importr('stats')
# poisson distribution
ppois = stats.ppois

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '$Revision: 442 $'

# Copyright 2010,2011 Jay R. Hesselberth

def identify_peaks(gdfilename, trackname, peak_width, chrom, strand,
                   log_pval_thresh, trim_peaks, shuffle_data,
                   local_lambda, max_width, verbose):

    warnings.simplefilter('ignore')
    with Genome(gdfilename) as genome:

        genome_coverage = calc_genome_coverage(genome, trackname)
      
        lambda_base = genome_coverage * peak_width

        if verbose:
            print >>sys.stderr, ">> Base lambda: %s" % lambda_base

        for chromosome in genome:

            # parallelize by chromosome if requested
            if chrom and chromosome.name != chrom: continue

            track_index = chromosome.index_continuous(trackname)

            # call peaks, merge them and trim if requested
            peak_data = identify_local_peaks(chromosome, track_index,
                                             peak_width, lambda_base,
                                             log_pval_thresh,
                                             local_lambda, 
                                             shuffle_data, verbose)

            peak_data = merge_adjacent_windows(peak_data, peak_width,
                                               verbose)

            if trim_peaks and not shuffle_data:
                peak_data = find_trimmed_peaks(peak_data, genome,
                                               chromosome,
                                               trackname, max_width,
                                               verbose)

            # output peaks for this chrom
            for start, end, pvalue in peak_data:
                fields = (chromosome.name, start, end,
                          trackname, pvalue, strand)
                print '\t'.join(map(str, fields))
                          
def merge_adjacent_windows(peak_data, peak_width, verbose):
    ''' merges adjacent significant windows that are peak_width apart. the
    more significant pvalue among the windows is returned.

    returns: merged (start, end, pval) tuples'''

    # TODO: recalculate pvalue for the merged window and return that

    if verbose:
        print >>sys.stderr, ">> merging windows ..."

    merged_windows = []
    last_start = last_end = last_pval = None
    adj_start = adj_end = adj_pval = None

    # sort data JIC, start to finish
    peak_data.sort()

    for start, end, pval in peak_data:

        # check to see whether current end and last_end are peak_width
        # apart. if so, set last to cur and continue
        if last_end and end - last_end == peak_width:

            # if the start of the adjacent windows is not set, then set it
            # to the start of the previous window
            if not adj_start: adj_start = last_start

            # update the end of the adjacent windows
            adj_end = end

            # set the pval for the adjacent widnows to the max pvalue btw
            # the current and last pvalues
            adj_pval = max(pval, last_pval, adj_pval)

        elif adj_start:
            # if the previous windows were adjacent, add them
            merged_windows.append((adj_start, adj_end, adj_pval))
            adj_start = adj_end = adj_pval = None

        elif last_start:
            # add the previous window
            merged_windows.append((last_start, last_end, last_pval))

        last_start, last_end, last_pval = start, end, pval

    # add the last remaining window
    if adj_start:
        merged_windows.append((adj_start, adj_end, adj_pval))
    else:
        merged_windows.append((last_start, last_end, last_pval))

    return merged_windows
    
def find_trimmed_peaks(peak_data, genome, chrom, trackname, max_width, verbose):
    ''' trim peaks to min width with max signal.

    returns: data with same form as peak_data, trimmed'''

    # TODO: recalculate pvalue for the trimmed and return that

    # XXX: there is likely a faster way to do this with native numpy
    # functions. speed scales exponentially with window size, may want to
    # add a "max_width" param to prevent trimming of peaks e.g.  >1000 bp

    trimmed_peaks = []

    if verbose:
        progress = ProgressBar(len(peak_data),
                               label=">> trimming peaks: ")

    for peak_start, peak_end, pval in peak_data:

        if peak_end - peak_start > max_width:
            trimmed_peaks.append((peak_start, peak_end, pval))
            continue

        if verbose: progress.next()

        continuous = genome[chrom.name][peak_start:peak_end, trackname]
        trim_peak_start, trim_peak_end = max_contiguous_signal(continuous)

        # update coords with new start & end - set end first
        peak_end = peak_start + trim_peak_end
        peak_start = peak_start + trim_peak_start

        trimmed_peaks.append((peak_start, peak_end, pval))

    # clean up progress bar
    if verbose: progress.end()

    return trimmed_peaks

def max_contiguous_signal(continuous):
    '''identify maximum contiguous signal from continuous data. 
    
    returns: start, end of contiguous region'''
    counts = []
    continuous_range = range(len(continuous)+1) # add 1 b/c half open

    max_count = 0
    # XXX: count_slop should be an option
    count_slop = 10
    # memoize nan positions
    nan_pos = dict()

    for start, end in combinations(continuous_range, r=2):

        # memoize nan positions - this should speed things up
        isnanpos = False
        for idx in (start, end-1):
            if idx in nan_pos:
                isnanpos = True
                break
            elif isnan(continuous[idx]):
                nan_pos[idx] = 1
                isnanpos = True

        if isnanpos: continue

        # make sure there is data in the peak
        score = nansum(continuous[start:end])

        if isnan(score): continue
        if score < max_count: continue

        counts.append((score, start, end))
       
        max_count = score

    # now find the smallest peak with the max count
    dists = [(end - start, start, end)
             for count, start, end in counts
             if count >= max_count - count_slop]

    dists.sort() # small to big
    trim_dist, trim_start, trim_end = dists[0]

    return (trim_start, trim_end)

def num_supercontigs(chromosome):
    ''' calculate number of Supercontigs in a genomedata.Chromosome '''
    num_supercontigs = 0
    for supercontig, _ in chromosome.itercontinuous():
        num_supercontigs += 1
    return num_supercontigs

def identify_local_peaks(chromosome, track_index, peak_width, lambda_base,
                         log_pval_thresh, local_lambda, shuffle_data, verbose):
    ''' identify peaks within a chromsome.
    
    returns: list of (start, end, pvalue) tuples'''

    local_peaks = []

    # pdb.set_trace()

    if verbose:
        progress = ProgressBar(num_supercontigs(chromosome),
                               label=">> calling peaks: ")

    for supercontig, continuous in chromosome.itercontinuous():

        if verbose: progress.next()

        # load data into memory
        track_continuous = continuous[:, track_index]

        # shuffle the data if requested
        if shuffle_data: shuffle(track_continuous)

        # generate window ranges i.e. peak widths
        peak_range = xrange(0, len(track_continuous), int(peak_width))

        for peak_start in peak_range:

            # set the peak end
            peak_end = int(peak_start + peak_width)

            peak_data = track_continuous[peak_start:peak_end]

            # skip the window if it is empty
            if isnan(nansum(peak_data)): continue

            # calculate local lambdas and choose max
            if local_lambda:
                lambda_values = calc_local_lambdas(peak_width,
                                                   track_continuous,
                                                   peak_start, peak_end)
                lambda_values.append(lambda_base)
                lambda_value = max(lambda_values)
            else:
                lambda_value = lambda_base

            # initial pvalue for the region
            log_pvalue = calc_log_pvalue(peak_data, lambda_value,
                                         verbose)

            # check that the pvalue is defined & that it meets the
            # threshold
            if not log_pvalue or log_pvalue < log_pval_thresh: continue

            # calculate chromosomal coords relative to current
            # contig coords
            chrom_start = supercontig.start + peak_start
            chrom_end = supercontig.start + peak_end

            # save current peak
            fields = (chrom_start, chrom_end, log_pvalue)
            local_peaks.append(fields)

            # reset lambda value to base
            lambda_value = lambda_base

    if verbose: progress.end() 

    return local_peaks

def calc_local_lambdas(peak_width, continuous, start, end):
    ''' calculate local lambdas similar to MACS "dynamic background.
    calculate lambada for different windows surround the candidate peak
    and return the local lambdas.
    
    returns: list of floats'''

    local_lambdas = []

    # XXX should be options
    local_widths = [1000, 5000, 10000]

    # calculate the width and center of the window
    window_width = end - start
    center = start + ((end - start) / 2)

    for local_width in local_widths:

        halfwidth = local_width / 2

        local_start = center - halfwidth
        local_end = center + halfwidth

        if local_start < 0 or local_end > len(continuous): continue

        # calculate the total signal and the number of defined positions
        total = nansum(continuous[local_start:local_end])

        lambda_local = total / local_width * peak_width 

        local_lambdas.append(lambda_local)

    return local_lambdas

def calc_log_pvalue(data, lambda_value, verbose):
    ''' calculate possion p-value for data in window, given a lambda value
    
    returns: abs(log10(pvalue)) (float) '''

    total = nansum(data)
    if isnan(total): return None

    # calculate p-value with rpy2.stats.ppois.  calculation is done in
    # log-space to avoid underflows
    # see:
    # http://blog.vladimirliu.com/2011/05/calculate-poisson-cdf-in-logarithm-space/

    pval = ppois(int(total), lambda_value, lower_tail=False, log_p=True)[0]

    return abs(pval)

def calc_genome_size(genome):
    ''' calculate genome size.

    returns: float'''
    genome_size = 0.0

    for chrom in genome:
        genome_size += chrom.end

    return genome_size

def calc_genome_coverage(genome, trackname):
    ''' calculate total number of reads mapping by summing signal and
    dividing by the number of positions with data
    
    returns: float'''

    total_coverage = 0.0

    for chrom in genome:
        track_idx = chrom.index_continuous(trackname)
        total_coverage += chrom.sums[track_idx]

    genome_size = calc_genome_size(genome)

    return total_coverage / genome_size

def parse_options(args):
    from optparse import OptionParser, OptionGroup
    
    description = ("Compute poisson p-values for peaks.")
    usage = '%prog [options] GENOMEDATADIR'
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
    group.add_option("-w", "--peak-width",
        action="store", type='float', dest="peak_width",
        help="peak width [default: %default]",
        metavar="WIDTH",
        default=150)

    group.add_option("-p", "--log-pvalue-threshold",
        action="store", type="float", dest="pval_thresh",
        help="log pvalue threshold [default: %default]",
        default=5) 

    group.add_option("-c", "--chrom",
        action="store", 
        help="parallelize by specifying chromosome name [default: %default]",
        default=None) 

    group.add_option("-s", "--strand",
        action="store", 
        help="strand in BED output [default: '%default']",
        default='.') 

    group.add_option("-M", "--max-width",
        action="store", 
        help="maximum width for window to be trimmed [default: '%default']",
        default=1000) 

    parser.add_option_group(group)

    group = OptionGroup(parser, "Flags")
    group.add_option("--trim-peaks",
        action="store_true", 
        help="trim initial peaks to their boundaries [default: %default]",
        default=False)    

    group.add_option("--local-lambda",
        action="store_true", 
        help="use local background for lambda [default: %default]",
        default=False)    

    group.add_option("--shuffle-data",
        action="store_true", 
        help="shuffle data for null p-values [default: %default]",
        default=False)    

    group.add_option("-v", "--verbose",
        action="store_true", dest="verbose",
        help="maximum verbosity [default: %default]",
        default=False)    

    parser.add_option_group(group)

    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error('specify GENOMEDATADIR')
    if not options.trackname:
        parser.error('specify trackname')

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    gdfilename = args[0]
    kwargs = {'trackname':options.trackname,
              'peak_width':options.peak_width,
              'log_pval_thresh':options.pval_thresh,
              'chrom':options.chrom,
              'strand':options.strand,
              'trim_peaks':options.trim_peaks,
              'max_width':options.max_width,
              'shuffle_data':options.shuffle_data,
              'local_lambda':options.local_lambda,
              'verbose':options.verbose}

    return identify_peaks(gdfilename, **kwargs)

if __name__ == '__main__':
    sys.exit(main())

