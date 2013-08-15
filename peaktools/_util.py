#! /usr/bin/envN python

'''peaktools._utils: common methods
'''

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '$Revision: 250 $'

# Copyright 2011, Jay R. Hesselberth

from itertools import product
from numpy import nansum

def enumerate_subpeaks(data, min_width=None):
    ''' return counts and coords of subpeaks contained within a parent
    peak.
    
    returns: tuples of (count, start, end)'''
    counts = []
    data_range = range(len(data))

    for start, end in product(data_range,repeat=2):
        if end < start: continue
        if min_width and abs(end - start) < min_width: continue
        counts.append((nansum(data[start:end]),start,end))

    return counts


