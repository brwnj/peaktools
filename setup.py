#!/usr/bin/env python

"""peaktools: tools for analyzing genomic peak calls"""

__version__ = "0.1"

# $Revision: 439 $

import os
import sys
import pdb
import pkg_resources

from ez_setup import use_setuptools
use_setuptools()

from platform import system, processor
from setuptools import setup

doclines = __doc__.splitlines()
name, short_description = doclines[0].split(": ")
long_description = "\n".join(doclines[2:])

# url = "XXXURL/%s/" % name.lower()
# download_url = "%s%s-%s.tar.gz" % (url, name, __version__)

classifiers = ["Natural Language :: English",
               "Programming Language :: Python"]

entry_points = """
[console_scripts]
peaktools-identify-peaks = peaktools.identify_peaks:main
peaktools-interolate-peaks = peaktools.interpolate_peaks:main
peaktools-qvalues = peaktools.qvalues:main
peaktools-qvalue-summary = peaktools.qvalue_summary:main
peaktools-color = peaktools.color:main
peaktools-combine-replicates = peaktools.combine_replicates:main
peaktools-test = test.run_tests:main
"""

install_requires = ["numpy", "rpy2", "textinput", "genomedata", "optplus",
                    "colorbrewer", "segtools"]

arch = "_".join([system(), processor()])

class InstallationError(Exception):
    pass

if __name__ == "__main__":
    setup(name=name,
          version=__version__,
          description=short_description,
          author="Jay Hesselberth",
          author_email="jay.hesselberth@gmail.com",
          # url=url,
          # download_url=download_url,
          classifiers=classifiers,
          long_description=long_description,
          install_requires=install_requires,
          zip_safe=False,
          # XXX: this should be based off of __file__ instead
          packages=['peaktools'], # 
          include_package_data = True,
          entry_points=entry_points
          )
