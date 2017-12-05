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


FILES STRUCTURE
===============
The files and folders included in Frigus are:

LICENSE                  licensing information of FRIGUS
MANIFEST.in              list of files that are included when generating the source distribution archive through setup.py
README.rst               this file
requirements.txt         list of packages that should be installed before using FRIGUS
setup.cfg                parameter configuration when building and installing FRIGUS through setup.py
setup.py                 script that is used to build and install FRIGUS
data/                    data directory needed to run examples and tests
src/                     source code of FRIGUS

The trees for src/ and data/ are described in the following, together with the files they contain.

-src/
It contains the folder python/, that is structured as follows:

├── examples/                                             folder containig the examples that can be run to calculate the cooling functions
│   ├── cooling_function_lipovka.py                       compute cooling function for the Lipovka cooling function
│   ├── cooling_function_three_level.py                   compute population densities for a 3-level system
│   ├── cooling_function_two_level.py                     compute population densities for a 2-level system
│   └── solve_equilibrium_time_dependent_paper_figure.py  time evolution example script of population densities
├── frigus                                                the "frigus" python package
│   ├── __init__.py                                       initialize module for the "frigus" python module
│   ├── metadata.py                                       metadata of "frigus", version number, authors...etc..
│   ├── population.py                                     calculate the level populations, cooling function
│   ├── read_collision_coefficients.py                    reader(s) for the collisional coefficients
│   ├── read_einstein_coefficient.py                      reader(s) for the Einstein coefficients
│   ├── read_energy_levels.py                             reader(s) for the energy levels of the multi-levels system
│   ├── readers.py                                        wrapper of the readers for the Einstein and collisional cofficients and energy levels
│   ├── utils.py                                          various utility functions
├── Makefile
└── tests
    ├── test_cooling_function_fits.py                     integration  tests for the cooling function fits and cooling function computations
    ├── test_solve_equilibium.py                          integration tests for a 2-level system
    └── test_two_level_system.py                          integration tests for a 3-level system

-data/
├── read/                                                 contains the data read from files
│   └── lipovka/                                          contains the data adopted for the cooling function calculation by Lipovka et al. 2005
│       ├── flower_roueff_data.dat                        HD rovibrational energies and collisional coefficients for the system H/HD
│       └── hd_einstein_coeffs.dat                        Einstein coefficients for HD as by Coppola et al. 2011
├── README                                                file describing the structure of the data/ folder
├── three_levels_1                                        contains the data to run FRIGUS in the case of a 3-levels system
│   ├── A_matrix.txt                                      the sample Einstein coefficients of the 3-levels system
│   ├── energy_levels.txt                                 the sample energy levels of the three-level system
│   ├── K_dex_matrix.txt                                  the sample collisional coeffients of the three-level system
│   └── README                                            information about the three-level system
└── two_levels_1:                                         contains the data to run Frigus in the case of a 2-levels system
    ├── A_matrix.txt                                      the sample Einstein coefficients of the 2-levels system
    ├── energy_levels.txt                                 the sample energy levels of the 2-levels system
    ├── K_dex_matrix.txt                                  the sample collisional coeffients of the 2-levels system
    └── README                                            information about the 2-levels system
