# resolvedvelocities
Python implementation of the Heinselman and Nicolls Bayesian reconstruction algorithm used to derive resolved ion drift velocity vectors

**Install**::

	git clone https://github.com/amisr/resolvedvelocities.git
	cd resolvedvelocities
	pip install .

**Run**:
Modify ``example_config.ini`` to include the correct file path to an input AMISR fitted file and the output file path and name you would like.
Then run::

	resolvedvelocities example_config.ini

**Note**:
This requires `apexpy <https://github.com/aburrell/apexpy>`_, which MUST be installed from commit 94e63c3af524b60dd652eb1c92df786e7ac9076a. This version includes the function ``bvectors_apex``, which returns the value of Be3, nessisary to convert the plasma drift velocity to the electric field. 

