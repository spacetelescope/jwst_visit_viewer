JWST Visit Viewer
--------------------------------------------------------------


This is a tool to display and visualize the pointings for JWST visit files,
such as are produced by OPGS based on APT files for upload to the observatory,
where they will be used by the OSS Executive to orchestrate observations.  The
main point of this tool is to make it easy to visualize where in the sky a
given visit will be pointed, the observatory attitude with respect to the sun
and the field of regard constraints, and which detector(s) will be used to take
data.

**Caveats, warnings, and disclaimers:** This is an unofficial tool provided on a best-effort basis.
Coordinate transforms may not yet be precisely correct in all details. Some aspects of V to J frame
transforms are not yet tracked in the SIAF PRD content.

A difference in the assumed ephemeris for JWST can result in inconsistencies in
field of regard calculations. If this tool warns that a given visit may be out
of the field of regard, take that with a grain of salt for now, and check with
the real PPS experts for a more authoritative answer.


Installation
------------


**Requirements:**
- numpy, astropy, matplotlib, etc
- pysiaf: https://github.com/spacetelescope/pysiaf
- jwst_gtvt: https://github.com/spacetelescope/jwst_gtvt

_Note, an up-to-date dev version of pysiaf is required, including this PR recently merged: https://github.com/spacetelescope/pysiaf/pull/177_


**Basic Installation**:

Install this repo:

    pip install git+https://github.com/spacetelescope/jwst_visit_viewer.git

Install dependencies (assuming you already have the basics):

    pip install git+https://github.com/spacetelescope/jwst_gtvt.git

    pip install git+https://github.com/spacetelescope/pysiaf.git


Usage Instructions
-------------------------------------------------

See [this notebook](https://github.com/spacetelescope/jwst_visit_viewer/blob/master/docs/Visit%20Viewer%20Docs%20and%20Usage.ipynb) in the Docs directory.


Contributing Code, Documentation, or Feedback
---------------------------------------------

See [CONTRIBUTING.md](CONTRIBUTING.md)

License
-------

BSD. See [LICENSE.md](LICENSE.md)
