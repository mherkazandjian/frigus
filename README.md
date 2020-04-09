[![Build Status](https://travis-ci.org/mherkazandjian/frigus.svg?branch=master)](https://travis-ci.org/mherkazandjian/frigus)

FRIGUS
======

Frigus is a package for computing the equilibrium and time-dependent level
populations of species. The level populations are used to compute the cooling
function that are parametrized in terms of the gas density and the kinetic and
radiation temperatures.

Obtaining
---------

The official source code repository is

   https://github.com/mherkazandjian/frigus

The source code can also be obtained from

   https://pypi.org/project/frigus/

or as wheels also from pypi.


Installation
------------

There are several way for installing frigus. For users the easiest way
is to install ``Frigus`` through ``pip`` via:

    pip install frigus

or using pip and directly from github

    pip install https://github.com/mherkazandjian/frigus/archive/master.zip --force --upgrade

Tests
-----

The tests can be run by executing (in the directory ``src/python``)

````bash
   make test
````

Input data
----------

The species data are located in the ``data`` directory. Each subdirectory
correspends to a certain set of data. The sub-directories ``data/two_levels_1``
and ``data/three_levels_1`` can be used as templates to construct datasets for
sophisticated species.


