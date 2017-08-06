Getting started
===============

Obtaining and Installation
--------------------------

The official repository of Frigus is:

    git@github.com:mherkazandjian/frigus.git

The latest version of Frigus can be installed through:

    pip install git+https://github.com/mherkazandjian/frigus.git@v3.0

The development version can be installed through

    pip install git+https://github.com/mherkazandjian/frigus.git@master

Alternatively the repository can be cloned and installed locally

    git clone https://github.com/mherkazandjian/frigus.git
    cd frigus
    python setup.py install

To run the tests after cloning the repo

    cd frigus/src/python/tests
    PYTHONPATH=..:${PYTHONPATH} py.test -v .

Building the documentation
--------------------------
The documentation is handled through the sphinx package and is
available at https://readthedocs.org/frigus ( .. todo:: set correct link)

    cd frigus/doc/sphinx
    PYTHONPATH=../../src/python:$PYTHONPATH make html
    firefox _build/html/index.html
