resolvedvelocities
==================

Overview
--------
`resolvedvelocities` is a Python implementation of the Heinselman and Nicolls Bayesian reconstruction algorithm [1]_ used to resolved 3D ion drift velocity and electric field vectors from AMISR line-of-sight measurements. In addition, this library uses a rotated pole transformation to avoid distorted latitude bins near the poles.

Quick Start
-----------
This package requires `numpy <https://numpy.readthedocs.io/en/latest/>`_ and `apexpy <https://apexpy.readthedocs.io/en/latest/>`_ and can be installed with `pip`::

	pip install git+https://github.com/amisr/resolvedvelocities.git

This code is designed to be run from the command line with a configuration file::

	resolvedvelocities config.ini

The configuration file specifies the input and output file names as well as a variety of other parameters.

**Note**:
The `apexpy <https://apexpy.readthedocs.io/en/latest/>`_ dependency MUST be installed from commit 94e63c3af524b60dd652eb1c92df786e7ac9076a. This version includes the function ``bvectors_apex``, which returns the value of Be3, necessary to convert the plasma drift velocity to the electric field.

Documentation
-------------

Full documentation is available at `<https://resolvedvelocities.readthedocs.io/en/latest/>`_.

Source Code
-----------

Source code is available from the `resolvedvelocities GitHub page <https://github.com/amisr/resolvedvelocities>`_.  Code is open source under GNU General Public License v3.0.  Please report bugs by submitting an `issue <https://github.com/amisr/resolvedvelocities/issues>`_!

References
----------

.. [1] Heinselman, C. J., and Nicolls, M. J. (2008), A Bayesian approach to electric field and E‐region neutral wind estimation with the Poker Flat Advanced Modular Incoherent Scatter Radar, Radio Sci., 43, RS5013, doi:`10.1029/2007RS003805 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007RS003805>`_.
