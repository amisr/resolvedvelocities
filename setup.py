#!/usr/bin/env python
"""
resolvedvelocities calculates the 3D plasma drift velocity
  and convection electric field from AMISR LoS velocities
  using the algorithm described by Heinselman and Nicolls, 2008
The full license can be found in LICENSE.txt
"""

import os
import re
import sys
import subprocess
from setuptools import find_packages, setup

here = os.path.abspath(os.path.dirname(__file__))

# Get the package requirements
REQSFILE = os.path.join(here, 'requirements.txt')
with open(REQSFILE, 'r') as f:
    REQUIREMENTS = f.readlines()
REQUIREMENTS = '\n'.join(REQUIREMENTS)

# Get the readme text
README = os.path.join(here, 'README.rst')
with open(README, 'r') as f:
    READMETXT = f.readlines()
READMETXT = '\n'.join(READMETXT)

# Get version number from __init__.py
regex = "(?<=__version__..\s)\S+"
with open(os.path.join(here,'resolvedvelocities/__init__.py'),'r', encoding='utf-8') as f:
    text = f.read()
match = re.findall(regex,text)
version = match[0].strip("'")

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
      version=version,
      description=DESC,
      author="AMISR",
      author_email="leslie.lamarche@sri.com",
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
              'resolvedalts=resolvedalts.ResolveVectors:main',
        ],
}
      )
