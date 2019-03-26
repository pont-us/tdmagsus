tdmagsus: tools for thermomagnetic susceptibility analysis
==========================================================

This is a small library for working with temperature-dependent magnetic
susceptibility data measured on an AGICO kappabridge (e.g. Agico, 2003;
Agico, 2009).

Overview
--------

tdmadsus provides three classes:

``Furnace`` represents the temperature-susceptibility behaviour of the empty
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
techniques (Petrovský  & Kapička, 2006), calculate the alteration index, and
write the data to a CSV file.

``MeasurementSet`` represents the data from a progressive sequence of
heating-cooling cycles, which it stores as a dictionary of ``MeasurementCycle``
objects indexed by peak temperature. It is initialized from a directory
containing multiple ``.CUR`` files.

License
-------

Copyright 2019 Pontus Lurcock; released under the `GNU General Public License,
version 3.0 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_. See the file
``COPYING`` for details.

References
----------

Agico, 2003. *KLY-3 / KLY-3S / CS-3 / CS-L / CS-23 user’s guide*, Brno,
Czech Republic: Advanced Geoscience Instruments Co.
https://www.agico.com/downloads/documents/manuals/kly3-man.pdf

Agico, 2009. *MFK1-FA / CS4 / CSL, MFK1-A / CS4 / CSL, MFK1-FB, MFK1-B user’s
guide* 4th ed., Brno, Czech Republic: Advanced Geoscience Instruments Co.
https://www.agico.com/downloads/documents/manuals/mfk1-man.pdf

Petrovský, E. & Kapička, A., 2006. On determination of the Curie point from
thermomagnetic curves. *Journal of Geophysical Research*, 111, p.B12S27.
