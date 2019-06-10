# resolvedvelocities
Python implementation of the Heinselman and Nicolls Bayesian reconstruction algorithm used to derive resolved ion drift velocity vectors

To Install::
	git clone https://github.com/amisr/resolvedvelocities.git
	cd resolvedvelocities
	pip install .

To Run:
Modify ``example_config.ini`` to include the correct file path to an input AMISR fitted file and the output file path and name you would like.
Then run::
	resolvedvelocities example_config.ini

Note:
This depends on apexpy, which MUST be installed from commit 94e63c3af524b60dd652eb1c92df786e7ac9076a.


