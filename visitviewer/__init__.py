# Visit Quick View: 
#   A toolkit for rapidly & intuitively understanding what a given OSS visit file will do


from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    __version__ = 'unknown'

import argparse
import os,sys
if not os.path.dirname(__file__) in sys.path:
    sys.path.append(os.path.dirname(__file__))
from . import visitparser, visitplotter
from .visitparser import VisitFileContents

def view(visitfilename, verbose=False, save=False, **kwargs):
    """top-level interface for visit quick view"""
    visit = visitparser.VisitFileContents(visitfilename, **kwargs)
    visitplotter.multi_plot(visit, save=save, verbose=verbose, **kwargs)

def main():
    """ Main function for command line arguments """
    parser = argparse.ArgumentParser(
        description='Visit Viewer'
    )
    parser.add_argument('filename', metavar='filename', type=str, help='Filename of visit file, or filenames for multiple', nargs='+')
    parser.add_argument('-d', '--dss', action='store_true', help='Use DSS instead of 2MASS')
    parser.add_argument('-n', '--no_gspa_yoffset', action='store_true', help='Do not apply FGS Yics angle offset to GSPA parameter. Use this for consistency with PPS version 14.14.1 and earlier (as used in LRE3, LRE4)')
    args = parser.parse_args()

    for filename in args.filename:
        view(filename, save=True, use_dss=args.dss, no_gspa_yoffset=args.no_gspa_yoffset)

if __name__ == '__main__':
    main()
