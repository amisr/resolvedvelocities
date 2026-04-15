resolvedvelocities
==================

Overview
--------
`resolvedvelocities` is a Python implementation of the Heinselman and Nicolls Bayesian reconstruction algorithm [1]_ used to resolved 3D ion velocity from AMISR line-of-sight measurements. The package contains two modules:

*iondrifts* calculates the ExB ion drift velocity and corresponding electric field at altitudes where the plasma is magnitized.  There are two flavors of this algorithm, one which bins points by magnetic latitude and is typically used in the F-region to get local plasma convection velocities, and one which bins point by altitude and is typically used in the E-region to get altitude profiles of velocity.

*neutralwinds* calculates the ion velocity, neutral wind velocity, and electric field in a altitude colum where the plasma is collisional.  This can typically only be used in the lower E-region, but requires data from higher altitude to estimate the electric field.

Quick Start
-----------
This package requires `numpy <https://numpy.readthedocs.io/en/latest/>`_ and `apexpy <https://apexpy.readthedocs.io/en/latest/>`_ and can be installed with `pip`::

	pip install git+https://github.com/amisr/resolvedvelocities.git

This code is designed to be run from the command line with a configuration file.  The configuration file specifies the input and output file names as well as a variety of other parameters.

Latitude-Binned Resolved Velocities (Vvels-Lat)::

	resolvedvelocities-lat config.ini

Altitude-Binned Resolved Velocities (Vvels-Alt)::

	resolvedvelocities-alt config.ini

Altitude-Binned Neutral Winds::

        resolvedvelocities-nw config.ini

Note that a different config file is currently needed for the neutral winds processing.


**Note**:
This package requires `numpy <https://numpy.readthedocs.io/en/latest/>`_ and `apexpy <https://apexpy.readthedocs.io/en/latest/>`_, which may not install correctly automatically on some systems.  If you encounter problems installing this package, try installing both of these packages independently first following their own installation instructions before attempting to reinstall `resolvedvelocities`.

Documentation
-------------

Full documentation is available at `<https://resolvedvelocities.readthedocs.io/en>`_.

Source Code
-----------

Source code is available from the `resolvedvelocities GitHub page <https://github.com/amisr/resolvedvelocities>`_.  Code is open source under GNU General Public License v3.0.  Please report bugs by submitting an `issue <https://github.com/amisr/resolvedvelocities/issues>`_!

References
----------

.. [1] Heinselman, C. J., and Nicolls, M. J. (2008), A Bayesian approach to electric field and E‐region neutral wind estimation with the Poker Flat Advanced Modular Incoherent Scatter Radar, Radio Sci., 43, RS5013, doi:`10.1029/2007RS003805 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007RS003805>`_.
