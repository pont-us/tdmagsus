tdmagsus: tools for thermomagnetic susceptibility analysis
==========================================================

This is a small library for working with temperature-dependent magnetic
susceptibility data measured on an AGICO kappabridge.

Overview
--------

tdmadsus provides three classes:

``Furnace``represents the temperature-susceptibility behaviour of the empty
furnace (i.e. the measurement apparatus without a sample). It allows a "raw"
set of sample measurements to be corrected to remove the effects of the
changes in the susceptibility of the equipment itself. A ``Furnace`` object is
created from a ``.CUR`` file produced from a measurement run with no sample.
Since measuremed data is frequently noisy, ``Furnace`` provides methods for
smoothing the data with a spline before it is used for corrections.

``MeasurementCycle`` represents the temperature-susceptibility behaviour of a
sample during a single heating-cooling sample. It is initialized from a ``.CUR``
file and, optionally, a ``Furnace`` object. If a furnace is supplied, it is
used to correct the measured sample data. ``MeasurementCycle`` provides methods
to calculate a disordering (Curie or Néel) temperature using different
techniques, calculate the alteration index, and write the data to a CSV file.

``MeasurementSet`` represents the data from a progressive sequence of
heating-cooling cycles, which it stores as a dictionary of ``MeasurementCycle``
objects indexed by peak temperature. It is initialized from a directory
containing multiple ``.CUR`` files.

License
-------

Copyright 2019 Pontus Lurcock; released under the `GNU General Public License,
version 3.0 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_. See the file
``COPYING`` for details.
