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
from matplotlib.backends.backend_pdf import PdfPages

def view(visitfilename, verbose=False, save=False, **kwargs):
    """top-level interface for visit quick view"""
    visit = visitparser.VisitFileContents(visitfilename, **kwargs)
    visitplotter.multi_plot(visit, save=save, verbose=verbose, **kwargs)

def multi_view(filenamelist, output_dir=None, **kwargs):
    filenamelist.sort()
    visit_ids = [os.path.splitext(os.path.basename(f))[0] for f in filenamelist]
    outname = f"{visit_ids[0]}_to_{visit_ids[-1]}_view.pdf"
    if output_dir is not None:
        outname = os.path.join(output_dir, outname)
    with PdfPages(outname) as pdf:
        for fn in filenamelist:
            view(fn, **kwargs)
            pdf.savefig()
    print(f"Output saved to {outname}")


def main():
    """ Main function for command line arguments """
    parser = argparse.ArgumentParser(
        description='Visit Viewer'
    )
    parser.add_argument('filename', metavar='filename', type=str, help='Filename of visit file, or filenames for multiple', nargs='+')
    parser.add_argument('-d', '--dss', action='store_true', help='Use DSS instead of 2MASS')
    parser.add_argument('-m', '--multipage', action='store_true', help='If running on multiple input files, generate one multi-page PDF instead of separate PDFs')
    parser.add_argument('-n', '--no_gspa_yoffset', action='store_true', help='Do not apply FGS Yics angle offset to GSPA parameter. Use this for consistency with PPS version 14.14.1 and earlier (as used in LRE3, LRE4)')
    parser.add_argument('-o', '--output_dir',  help='Output PDF to this directory (rather than current working dir)')
    args = parser.parse_args()

    if args.multipage and len(args.filename) > 1:
        # Generate one multi page PDF output
        multi_view(args.filename, use_dss=args.dss, no_gspa_yoffset=args.no_gspa_yoffset, output_dir=args.output_dir)

    else:
        # Generate a separate PDF per each input
        for filename in args.filename:
            view(filename, save=True, use_dss=args.dss, no_gspa_yoffset=args.no_gspa_yoffset, output_dir=args.output_dir)

if __name__ == '__main__':
    main()
