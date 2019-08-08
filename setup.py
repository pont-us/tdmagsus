#!/usr/bin/env python3

import os.path
from setuptools import setup

cwd = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(cwd, "README.rst")) as fh:
    long_desc = fh.read()

setup(name="tdmagsus",
      version="1.0.2",
      description=
      "Manipulation of temperature-dependent magnetic susceptibility data",
      long_description_content_type="text/x-rst",
      long_description=long_desc,
      url="https://github.com/pont-us/tdmagsus",
      author="Pontus Lurcock",
      author_email="pont@talvi.net",
      license="GNU GPLv3+",
      classifiers=["Development Status :: 5 - Production/Stable",
                   "License :: OSI Approved :: "
                   "GNU General Public License v3 or later (GPLv3+)",
                   "Topic :: Scientific/Engineering",
                   "Programming Language :: Python :: 3",
                   "Intended Audience :: Science/Research"
                   ],
      packages=["tdmagsus"],
      install_requires=["numpy", "scipy"],
      zip_safe=False)
