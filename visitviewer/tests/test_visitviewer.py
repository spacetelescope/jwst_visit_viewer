import os.path
import glob

import visitviewer

EXAMPLES_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'examples')

def test_parser():
    """ Test we can at least parse all the example files
    """
    example_files = glob.glob(os.path.join(EXAMPLES_PATH, "V*vst"))

    for filename in example_files:
        print(filename)
        result = visitviewer.VisitFileContents(filename)


def test_plotter():
    """ Test we can make at least one plot
    """
    example_files = glob.glob(os.path.join(EXAMPLES_PATH, "V*vst"))

    result = visitviewer.view(example_files[0])


