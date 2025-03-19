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


def test_plot_field_of_regard_date():
    import os
    result = visitviewer.visitplotter.plot_field_of_regard_on_date('2024-12-25')
    output_fn_name = 'jwst_field_of_regard_2024-12-25.png'
    assert os.path.exists(output_fn_name)
    os.remove(output_fn_name)

