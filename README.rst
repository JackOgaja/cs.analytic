=========
Overview
=========

HOS Utilities

Analytic Utilities for developing and analysing Higher Order Conservative Discretization Methods 
used in Geophysical Modeling and Fluid Mechanics studies.

This documentation is under development.

Installation:
=============

To install, type the following command::

  ~$ git clone https://github.com/JackOgaja/cs.analytic.git 
  ~$ cd cs.analytic

Uninstallation:
===============

To uninstall, type the following command::

  ~$ rm -rf cs.analytic

Run and profile:
================
To run and profile the scripts type the following command::

  ~$ profile tests/tests.py

OR::

  ~$ python -m cProfile tests/tests.py

OR, import cProfile for individual functipons::

  ~$ import cProfile
  ~$ cProfile.run('main()')

Features:
=========

The following analysis can be performed:
   #. Phase and group velocity errors analysis
   #. Linear stability analysis
   #. Aliasing (Non-linear stability) analysis

License:
========

MIT License

For detailed copyright and licensing information, please refer to the
license file `LICENSE.md` in the top level cs.analytic directory.

