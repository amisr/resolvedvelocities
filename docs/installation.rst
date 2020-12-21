Installation
============

Basic Installation
------------------
`resolvedvelocites` is a pure python package so installation with `pip` is straightforwards::

  pip install git+https://github.com/amisr/resolvedvelocities.git

This code depends on `apexpy <https://apexpy.readthedocs.io/en/latest/>`_, which may have to be manually installed BEFORE you can install `resolvedvelocies`.  You MUST install a version of apexpy that includes commit 94e63c3af524b60dd652eb1c92df786e7ac9076a. This commit adds the function ``bvectors_apex``, which returns the value of Be3, necessary to convert the plasma drift velocity to the electric field.  To install from this commit, the standard pip install command can be modified as follows::
  
  pip install git+https://github.com/aburrell/apexpy.git@94e63c3af524b60dd652eb1c92df786e7ac9076a


Developers Installation
-----------------------

If you would like to contribute to or expand `resolvedvelocities`, clone the source repository and install from your new local copy::

  git clone https://github.com/amisr/resolvedvelocities.git
	cd resolvedvelocities
	pip install .
