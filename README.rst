[![Build Status](https://travis-ci.org/mherkazandjian/frigus.svg?branch=master)](https://travis-ci.org/mherkazandjian/frigus)

FRIGUS
======

Package for computing the equilibrium and time-dependent level populations of
species. The level populations are used to compute the cooling function that
are parametrized in terms of the gas density and the kinetic and radiation
temperatures.

Obtaining
---------

https://github.com/mherkazandjian/frigus

Installation
------------

There are several way for installing frigus. For users the easiest way
is to install ``Frigus`` through ``pip`` via:

    pip install frigus


Tests
-----

The tests can be run by executing

   make test

in the directory ``src/python``.

Input data
----------

The species data are located in the ``data`` directory. Each subdirectory
correspends to a certain set of data. The sub-directories ``data/two_levels_1``
and ``data/three_levels_1`` can be used as tempates to construct datasets for
sophisticated species.


