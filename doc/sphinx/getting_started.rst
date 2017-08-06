Getting started
===============

Obtaining and Installation
--------------------------

The official repository of Frigus is:

.. code-block:: bash

    git@github.com:mherkazandjian/frigus.git

The latest version of Frigus can be installed through:

.. code-block:: bash

    pip install git+https://github.com/mherkazandjian/frigus.git@v3.0

The development version can be installed through

.. code-block:: bash

    pip install git+https://github.com/mherkazandjian/frigus.git@master

Alternatively the repository can be cloned and installed locally

.. code-block:: bash

    git clone https://github.com/mherkazandjian/frigus.git
    cd frigus
    python setup.py install

To run the tests after cloning the repo

.. code-block:: bash

    cd frigus/src/python/tests
    PYTHONPATH=..:${PYTHONPATH} py.test -v .

Once Frigus is installed or the path to Frigus is set in the PYTHONPATH,
it can be imported through

.. code-block:: python

    import frigus


Building the documentation
--------------------------
The documentation is handled through the sphinx package and is
available at https://readthedocs.org/frigus ( .. todo:: set correct link)

.. code-block:: bash

    cd frigus/doc/sphinx
    PYTHONPATH=../../src/python:$PYTHONPATH make html
    firefox _build/html/index.html
