#! /usr/bin/env python

'''peaktools-color: select color for track'''

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '$Revision: 294 $'

# Copyright 2010,2011 Jay R. Hesselberth

import sys
import colorbrewer

# XXX mmh: why make colorbrewer schemes accessbile through globals?  why
# not a dict instead? much more flexible accessor
NUM_COLORS = 8
SCHEME = colorbrewer.Dark2[NUM_COLORS]

def _choose_color(discrete):
    color = SCHEME[discrete % NUM_COLORS]
    return ','.join(map(str, color))

def print_color(discrete):
    print _choose_color(discrete)

def parse_options(args):
    from optparse import OptionParser
    
    description = "Report RGB value string for discrete number"
    usage = '%prog [options] NUMBER'
    version = '%%prog %s' % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    options, args = parser.parse_args()

    if len(args) < 1:
        parser.error('specify discrete number')

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    discrete = int(args[0])
    kwargs = {}
    return print_color(discrete, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
  
