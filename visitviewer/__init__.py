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
    visit = visitparser.VisitFileContents(visitfilename)
    visitplotter.multi_plot(visit, save=save, verbose=verbose, **kwargs)

def main():
    """ Main function for command line arguments """
    parser = argparse.ArgumentParser(
        description='Visit Viewer'
    )
    parser.add_argument('filename', metavar='filename', type=str, help='Filename of visit file')
    args = parser.parse_args()
    view(args.filename, save=True)

if __name__ == '__main__':
    main()
