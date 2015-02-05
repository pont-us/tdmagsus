#!/usr/bin/env python3

from setuptools import setup

setup(name = "tdmagsus",
      version = "0.0",
      description =
      "Manipulation of temperature-dependent "
      "magnetic susceptibility data",
      url = "https://github.com/pontu/tdmagsus",
      author = "Pontus Lurcock",
      author_email = "pont@talvi.net",
      license = "GNU GPLv3+",
      classifiers = ["Development Status :: 4 - Beta",
                     "License :: OSI Approved :: "
                     "GNU General Public License v3 or later (GPLv3+)",
                     "Topic :: Scientific/Engineering",
                     "Programming Language :: Python",
                     "Intended Audience :: Science/Research`"
                 ],
      packages = ["tdmagsus"],
      install_requires = ["numpy", "scipy"],
      zip_safe = False)
