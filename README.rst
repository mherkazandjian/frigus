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
LICENSE
MANIFEST.in
README.rst
requirements.txt
setup.cfg
setup.py
data/
src/

The tree for src/ and data/ are described in the following, together with the files they contain.

-src/
It contains the folder python/, that is structured as follows:
├── examples: it contains the examples that can be run to calculate the cooling function; in particular they are:
│   ├── cooling_function_lipovka.py
│   ├── cooling_function_three_level.py
│   ├── cooling_function_two_level.py
│   └── solve_equilibrium_time_dependent_paper_figure.py
├── frigus: it contains the source codes:
│   ├── __init__.py
│   ├── __init__.pyc
│   ├── metadata.py
│   ├── metadata.pyc
│   ├── population.py: it calculates the level population of a multi-levels system in the steady state approximation
│   ├── population.pyc
│   ├── __pycache__
│   │   ├── __init__.cpython-36.pyc
│   │   ├── population.cpython-36.pyc
│   │   ├── read_collision_coefficients.cpython-36.pyc
│   │   ├── read_einstein_coefficient.cpython-36.pyc
│   │   ├── read_energy_levels.cpython-36.pyc
│   │   ├── readers.cpython-36.pyc
│   │   └── utils.cpython-36.pyc
│   ├── read_collision_coefficients.py: it reads the collisional coefficients introduced in the kinetic scheme
│   ├── read_collision_coefficients.pyc
│   ├── read_einstein_coefficient.py: it reads the Einstein coefficients introduced in the kinetic scheme
│   ├── read_einstein_coefficient.pyc
│   ├── read_energy_levels.py: it reads the energy levels of the multi-levels system
│   ├── read_energy_levels.pyc
│   ├── readers.py: it implements the readers for the provided examples
│   ├── readers.pyc
│   ├── utils.py
│   └── utils.pyc
├── Makefile
└── tests
    ├── htmlcov
    │   ├── coverage_html.js
    │   ├── _home_carla_science_mher_cooling_function_src_python_frigus___init___py.html
    │   ├── _home_carla_science_mher_cooling_function_src_python_frigus_population_py.html
    │   ├── _home_carla_science_mher_cooling_function_src_python_frigus_read_collision_coefficients_py.html
    │   ├── _home_carla_science_mher_cooling_function_src_python_frigus_read_einstein_coefficient_py.html
    │   ├── _home_carla_science_mher_cooling_function_src_python_frigus_read_einstien_coefficient_py.html
    │   ├── _home_carla_science_mher_cooling_function_src_python_frigus_read_energy_levels_py.html
    │   ├── _home_carla_science_mher_cooling_function_src_python_frigus_readers_py.html
    │   ├── _home_carla_science_mher_cooling_function_src_python_frigus_utils_py.html
    │   ├── index.html
    │   ├── jquery.ba-throttle-debounce.min.js
    │   ├── jquery.hotkeys.js
    │   ├── jquery.isonscreen.js
    │   ├── jquery.min.js
    │   ├── jquery.tablesorter.min.js
    │   ├── keybd_closed.png
    │   ├── keybd_open.png
    │   ├── status.json
    │   ├── style.css
    │   ├── test_cooling_function_fits_py.html
    │   ├── test_cooling_functions_py.html
    │   ├── test_fortran_py.html
    │   ├── test_solve_equilibium_py.html
    │   └── test_two_level_system_py.html
    ├── __pycache__
    │   ├── test_cooling_function_fits.cpython-27-PYTEST.pyc
    │   ├── test_cooling_function_fits.cpython-36-PYTEST.pyc
    │   ├── test_cooling_functions.cpython-27-PYTEST.pyc
    │   ├── test_fortran.cpython-27-PYTEST.pyc
    │   ├── test_solve_equilibium.cpython-27-PYTEST.pyc
    │   ├── test_solve_equilibium.cpython-36-PYTEST.pyc
    │   ├── test_two_level_system.cpython-27-PYTEST.pyc
    │   └── test_two_level_system.cpython-36-PYTEST.pyc
    ├── test_cooling_function_fits.py
    ├── test_solve_equilibium.py
    └── test_two_level_system.py
 

-data/
├── read: it contains the data read from files
│   └── lipovka: it contains the data adopted for the cooling function calculation by Lipovka et al. 2005
│       ├── flower_roueff_data.dat: HD rovibrational energies and collisional coefficients for the system H/HD
│       └── hd_einstein_coeffs.dat: Einstein coefficients for HD as by Coppola et al. 2011
├── README: file describing the structure of the data/ folder
├── three_levels_1: it contains the data to run Frigus in the case of a 3-levels system
│   ├── A_matrix.txt
│   ├── energy_levels.txt
│   ├── K_dex_matrix.txt
│   └── README
└── two_levels_1: it contains the data to run Frigus in the case of a 2-levels system
    ├── A_matrix.txt
    ├── energy_levels.txt
    ├── K_dex_matrix.txt
    └── README
 

 


