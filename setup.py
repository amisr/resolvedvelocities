#!/usr/bin/env python
""" 
resolvedvelocities calculates the 2D plasma drift velocity
  and convection electric field from AMISR LoS velocities
  using the algorithm described by Heinselman and Nicolls, 2008
The full license can be found in LICENSE.txt
"""

import os
import sys
import subprocess
from setuptools import find_packages, setup

# Get the package requirements
REQSFILE = os.path.join(os.path.dirname(__file__), 'requirements.txt')
with open(REQSFILE, 'r') as f:
    REQUIREMENTS = f.readlines()
REQUIREMENTS = '\n'.join(REQUIREMENTS)

# Do some nice things to help users install on conda.
if sys.version_info[:2] < (3, 0):
    EXCEPTION = OSError
else:
    EXCEPTION = subprocess.builtins.FileNotFoundError
try:
    subprocess.call(['conda', 'install', ' '.join(REQUIREMENTS)])
    REQUIREMENTS = []
except EXCEPTION:
    pass

# Get the readme text
README = os.path.join(os.path.dirname(__file__), 'README.rst')
with open(README, 'r') as f:
    READMETXT = f.readlines()
READMETXT = '\n'.join(READMETXT)

# Package description
DESC = "Tool for resolving 2D plasma drift from AMISR LoS velocity "

#############################################################################
# First, check to make sure we are executing
# 'python setup.py install' from the same directory
# as setup.py (root directory)
#############################################################################
PATH = os.getcwd()
assert('setup.py' in os.listdir(PATH)), \
       "You must execute 'python setup.py install' from within the \
repo root directory."


#############################################################################
# Now execute the setup
#############################################################################
setup(name='resolvedvelocities',
      install_requires=REQUIREMENTS,
      setup_requires=REQUIREMENTS,
      version="2.0.0",
      description=DESC,
      author="Ashton S. Reimer",
      author_email="ashtonsethreimer@gmail.com",
      url="https://github.com/amisr/resolvedvelocities",
      download_url="https://github.com/amisr/resolvedvelocities",
      packages=find_packages(),
      long_description=READMETXT,
      zip_safe=False,
      py_modules=['resolvedvelocities'],
      classifiers=["Development Status :: 2.0.0 - Release",
                   "Topic :: Scientific/Engineering",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Natural Language :: English",
                   "Programming Language :: Python",
                  ],
      entry_points={
          'console_scripts': [
              'resolvedvelocities=resolvedvelocities.run_resolve_vectors:main',
        ],
}
      )
