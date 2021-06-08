JWST Visit Viewer
--------------------------------------------------------------


This is a tool to display and visualize the pointings for JWST visit files,
such as are produced by OPGS based on APT files for upload to the observatory,
where they will be used by the OSS Executive to orchestrate observations.  The
main point of this tool is to make it easy to visualize where in the sky a
given visit will be pointed, the observatory attitude with respect to the sun
and the field of regard constraints, and which detector(s) will be used to take
data.

**Caveats, warnings, and disclaimers:** This is an official tool provided on a best-effort basis.
Coordinate transforms may not yet be precisely correct in all details. Some aspects of V to J frame
transforms are not yet tracked in the SIAF PRD content.


Installation
------------


**Requirements:**
- numpy, astropy, matplotlib, etc
- pysiaf: https://github.com/spacetelescope/pysiaf
- jwst_gtvt: https://github.com/spacetelescope/jwst_gtvt

_Note, currently a dev version of pysiaf is required, including this PR: https://github.com/spacetelescope/pysiaf/pull/177_


**Basic Installation**:

Install this repo:

    pip install git+https://github.com/spacetelescope/jwst_visit_viewer.git

Install dependencies (assuming you already have the basics):

    pip install git+https://github.com/spacetelescope/jwst_gtvt.git

    # Special, temporary: install dev branch with sky transforms:
    pip install git+https://github.com/mperrin/pysiaf.git@sky_transforms


Usage Instructions
-------------------------------------------------

See [this notebook](https://github.com/spacetelescope/jwst_visit_viewer/blob/master/docs/Visit%20Viewer%20Docs%20and%20Usage.ipynb) in the Docs directory.


Contributing Code, Documentation, or Feedback
---------------------------------------------

See [CONTRIBUTING.md](CONTRIBUTING.md)

License
-------

BSD. See [LICENSE.md](LICENSE.md)
