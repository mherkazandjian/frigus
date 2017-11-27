FRIGUS
======

Frigus is a package for computing the level populations of species.
The level populations can be used to compute cooling/heating functions.
In the current implementation, the following features are available:

   - equilibrium level populations
   - time-dependent evolution
   - cooling function

Obtaining
---------

Frigus can be obtained from:

   https://github.com/mherkazandjian/frigus.git

Installation
------------

   $ git clone https://github.com/mherkazandjian/frigus.git

   $ pip install -r requirements.txt

   $ python setup.py install

Tests
-----

The test suite can be run through:

   $ FRIGUS_SOURCE_ROOT/src/python>  make test

To run the examples:

   $ FRIGUS_SOURCE_ROOT/src/python/examples>  python cooling_function_two_level.py

   $ FRIGUS_SOURCE_ROOT/src/python/examples>  python cooling_function_three_level.py

The following two scripts repreduce two of the figures of the FRIGUS paper:

   $ FRIGUS_SOURCE_ROOT/src/python/examples>  python cooling_function_lipovka.py

   $ FRIGUS_SOURCE_ROOT/src/python/examples>  python solve_equilibrium_time_dependent_paper_figure.py




