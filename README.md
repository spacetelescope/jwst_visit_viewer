Space Telescope Science Institute Open Source Package Template
--------------------------------------------------------------
(put information here about the goal of this package and who uses it)

This is a standard package template for repositories under the "spacetelescope" organization. It includes a SPHINX documentation template that applies the STSCI branding to documentation and can be used to build documentation on ReadTheDocs or locally. 

This package uses [setuptools_scm](https://github.com/pypa/setuptools_scm) to create the version number. In order for the full build to show the correct version number in the documentation, assign a version tag to your repo, then install the package using `pip install .` or `python setup.py install`. Now you can go to the docs directory and build the documentation using "make html".

**If you are making heavy use of astropy, and would like to use astropy-helpers, use the [stsci-cookiecutter](https://github.com/spacetelescope/stsci-package-template/tree/stsci-cookiecutter) branch on this repository instead.**


Status reports for developers
-----------------------------
These are example tags that show the status of this repository for testing and documentation purposes
You can replace them with the badges for your package when testing, docs, and coverage are up and running.

> .. image:: https://travis-ci.org/spacetelescope/stsci-package-template.svg
>     :target: https://travis-ci.org/spacetelescope/stsci-package-template
>     :alt: Travis Status
>
> .. image:: https://readthedocs.org/projects/stsci-package-template/badge/?version=latest
>     :target: https://readthedocs.org/projects/stsci-package-template/?badge=latest
>     :alt: Documentation Status
>
> .. image:: https://coveralls.io/repos/github/spacetelescope/stsci-package-template/badge.svg?branch=master
>     :target: https://coveralls.io/github/spacetelescope/stsci-package-template?branch=master
>     :alt: Test Coverage Status


Installation
------------
Give an example of how to install this package


Where to find examples of how to use this package
-------------------------------------------------
If you choose to build documentation on [Readthedocs](https://readthedocs.org/), you should update the .rtd-environment.yml file to include all of your package dependencies. Some basic ones are already included. 


Testing Information
-------------------
The included .travis.yml file is a basic structure to run testing on [Travis-CI](https://docs.travis-ci.com/), replace occurances of "packagename" with the name of your software and add basic tests to the packagename/tests/ directory. There is a very basic example test in place already. 


Contributing Code, Documentation, or Feedback
---------------------------------------------



3rd Party Libraries this package requires
-----------------------------------------



License
-------
