#!/usr/bin/env python3

from setuptools import setup

setup(name="tdmagsus",
      version="1.0",
      description=
      "Manipulation of temperature-dependent magnetic susceptibility data",
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
