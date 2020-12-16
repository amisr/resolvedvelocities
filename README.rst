resolvedvelocities
==================

Overview
--------
`resolvedvelocities` is the python implementation of the Heinselman and Nicolls Bayesian reconstruction algorithm [1]_ used to resolved 3D ion drift velocity and electric field vectors from AMISR line-of-sight measurements.

Quick Start
-----------
This package requires `numpy <https://numpy.readthedocs.io/en/latest/>`_ and `apexpy <https://apexpy.readthedocs.io/en/latest/>`_ and can be installed with `pip`::

	pip install git+https://github.com/amisr/resolvedvelocities.git

This code is designed to be run from the command line with a configuration file::

	resolvedvelocities config.ini

The configuration file specifies the input and output file names as well as a variety of other parameters.

Documentation
-------------

Full documentation is available at <readthedocs link>

References
----------

.. [1] Heinselman, C. J., and Nicolls, M. J. (2008), A Bayesian approach to electric field and E‐region neutral wind estimation with the Poker Flat Advanced Modular Incoherent Scatter Radar, Radio Sci., 43, RS5013, doi:`10.1029/2007RS003805 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007RS003805>`_.

**Install**::

	git clone https://github.com/amisr/resolvedvelocities.git
	cd resolvedvelocities
	pip install .

**Run**:
Modify ``example_config.ini`` to include the correct file path to an input AMISR fitted file and the output file path and name you would like.
Then run::

	resolvedvelocities example_config.ini

**Note**:
This requires `apexpy <https://github.com/aburrell/apexpy>`_, which MUST be installed from commit 94e63c3af524b60dd652eb1c92df786e7ac9076a. This version includes the function ``bvectors_apex``, which returns the value of Be3, necessary to convert the plasma drift velocity to the electric field.
