#!/usr/bin/python3

# Copyright 2019 Pontus Lurcock.
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Tools for working with temperature-dependent magnetic susceptibility data.

This module contains classes for processing temperature-dependent magnetic
susceptibility data from AGICO kappabridges. It reads the CUR data
files produced by the kappabridge software and provides the ability to
subtract an empty furnace measurement (smoothed with a spline) from
sample measurements in order to provide a corrected measurement of the
susceptibility of the sample itself.

"""

import glob
import os.path
import re

import numpy
from numpy import array
from scipy.interpolate import UnivariateSpline

line_pattern = re.compile(r'^ +\d')
field_separator = re.compile(' +')


def read_cur_file(filename):
    """Read a .CUR magnetic susceptibility file.

    Args:
      filename (str): name of file to read

    Returns:
      ((heating_temps, heating_mag_sus_values],
       (cooling_temps, cooling_mag_sus_values])

      This is a tuple of two tuples, each containing
      two numpy arrays.
    """

    infile = open(filename, 'r')
    heating = ([], [])
    cooling = ([], [])
    current = heating
    cool = False
    prev_temp = -300
    for line in infile:
        if line_pattern.match(line.rstrip()):
            temperature, mag_sus = \
                map(float, field_separator.split(line.lstrip())[0:2])
            if temperature < prev_temp - 0.5 and not cool:
                current = cooling
                cool = True
            current[0].append(temperature)
            current[1].append(mag_sus)
            prev_temp = temperature
    infile.close()
    cooling[0].reverse()
    cooling[1].reverse()
    heating = (heating[0][1:], heating[1][1:],)
    cooling = (cooling[0][1:], cooling[1][1:],)
    return tuple(map(array, heating)), tuple(map(array, cooling))


class Furnace:

    @staticmethod
    def extend_data(temps_mss):
        """Pad data at ends to improve spline fit.

        Given ([T1, T2, T3, ... , Tn-1, Tn], [M1, M2, M3, ... , Mn-1, Mn]),
        produces ([T1-20, T1-10, T1, T2, T3, ... , Tn-1, Tn, Tn+10, Tn+20],
                  [M1, M1, M1, M2, M3, ... , Mn-1, Mn, Mn, Mn]).

        Args:
          temps_mss: tuple of temperatures and mag sus values

        Returns:
          same values, padded by two extra data points at each end
        """
        temps, mss = temps_mss
        tlist = temps.tolist()
        mlist = mss.tolist()
        tlist = [(tlist[0] - 20)] + [(tlist[0] - 10)] + \
            tlist + [tlist[-1] + 10] + [tlist[-1] + 20]
        mlist = [mlist[0]] + [mlist[0]] + \
            mlist + [mlist[-1]] + [mlist[-1]]
        return array(tlist), array(mlist)

    def __init__(self, filename, smoothing=100):
        """Initialize Furnace object from a CUR file.

        Args:
          filename: file path from which to read furnace data
          smoothing: smoothing factor for spline curve
        """
        heat, cool = read_cur_file(filename)
        self.heat_data, self.cool_data = map(Furnace.extend_data, (heat, cool))
        self.heat_spline = UnivariateSpline(*self.heat_data, s=smoothing)
        self.cool_spline = UnivariateSpline(*self.cool_data, s=smoothing)

    def get_spline_data(self):
        """Return furnace temperature/M.S. data and spline approximations.

        This method is mainly intended for checking that the splines are
        doing a good job of smoothing the data.

        Returns:
           (heating_data, heating_spline, cooling_data, cooling_spline)
           Each element of this tuple is itself a 2-tuple containing a 
           list of temperatures and a list of associated M.S. values."""

        splinex = numpy.arange(20, 701)
        spliney_heat = self.heat_spline(splinex)
        spliney_cool = self.cool_spline(splinex)
        return self.heat_data, (splinex, spliney_heat),\
            self.cool_data, (splinex, spliney_cool)    

    @staticmethod
    def correct_with_spline(temps, mss, spline):
        mss_corrected = numpy.zeros_like(mss)
        for i in range(0, len(temps)):
            mss_corrected[i] = mss[i] - spline(temps[i])
        return temps, mss_corrected

    def correct(self, heating, cooling):
        return (Furnace.correct_with_spline(heating[0], heating[1],
                                            self.heat_spline),
                Furnace.correct_with_spline(cooling[0], cooling[1],
                                            self.cool_spline))


class MeasurementCycle:
    """The results of a single heating-cooling run."""

    def __init__(self, furnace, filename, real_vol, nom_vol):
        self.furnace = furnace
        self.real_vol = real_vol
        self.nom_vol = nom_vol
        (heating, cooling) = read_cur_file(filename)
        if self.furnace is not None:
            heating, cooling = self.furnace.correct(heating, cooling)
        # heating = (heating[0], TdmsData.shunt_up(heating[1]))
        # cooling = (cooling[0], TdmsData.shunt_up(cooling[1]))
        heating = (heating[0], self.correct_for_volume(heating[1]))
        cooling = (cooling[0], self.correct_for_volume(cooling[1]))
        self.data = (heating, cooling)

    def write_csv(self, filename):
        """Write furnace-corrected data to a CSV file.

        Args:
          filename (str): name of file to write to.
        """

        cooling = (list(reversed(self.data[1][0])),
                   list(reversed(self.data[1][1])))
        with open(filename, "w") as fh:
            for direction in self.data[0], cooling:
                pairs = zip(direction[0], direction[1])
                for pair in pairs:
                    fh.write("%.2f,%.2f\n" % (pair[0], pair[1]))

    @staticmethod
    def chop_data(temps_mss, min_temp, max_temp):
        """Truncate data to a given temperature range.

        """
        temps, mss = temps_mss
        temps_out = []
        mss_out = []
        for i in range(0, len(temps)):
            temp = temps[i]
            if min_temp <= temp <= max_temp:
                temps_out.append(temp)
                mss_out.append(mss[i])
        return array(temps_out), array(mss_out)

    @staticmethod
    def linear_fit(xs, ys):
        fit = numpy.polyfit(xs, ys, 1)
        poly = numpy.poly1d(fit.tolist())
        model_ys = poly(xs)
        mean_y = numpy.mean(ys)
        sserr = numpy.sum((ys - model_ys)**2)
        sstot = numpy.sum((ys - mean_y)**2)
        rsquared = 1 - sserr / sstot
        return poly, rsquared

    def curie_paramag(self, min_temp, max_temp):
        """Estimate Curie temperature by linear fit to inverse susceptibility.

        Args:
          min_temp: minimum of temperature range for fit
          max_temp: maximum of temperature range for fit

        Return:
          (curie, rsquared, poly)
          curie: estimated Curie temperature
          rsquared: R² value for fit
          poly: polynomial object representing line of best fit
        """

        temps, mss = \
            MeasurementCycle.chop_data(self.data[0], min_temp, max_temp)
        poly, rsquared = MeasurementCycle.linear_fit(temps, 1./mss)
        curie = poly.r[0]  # x axis intercept
        return curie, rsquared, poly

    def curie_inflection(self, min_temp, max_temp):
        """Estimate Curie temperature by inflection point.

        Estimate Curie point by determining the inflection point of 
        the curve segment starting at the Hopkinson peak. The curve
        segment must be specified.

        Args:
          min_temp: start of curve segment
          max_temp: end of curve segment

        Result:
          (temp, spline)
          temp is the estimated Curie temperature
          spline is the scipy.interpolate.UnivariateSpline used to fit
            the data and determine the inflection point
        """

        # Fit a cubic spline to the data. Using the whole dataset gives
        # a better approximation at the endpoints of the selected range.
        spline = UnivariateSpline(self.data[0][0], self.data[0][1], s=.1)

        # Get the data points which lie within the selected range.
        temps, _ = MeasurementCycle.chop_data(self.data[0], min_temp, max_temp)

        # Evaluate the second derivative of the spline at each selected
        # temperature step.
        derivs = [spline.derivatives(t)[2] for t in temps]

        # Fit a new spline to the derivatives in order to calculate the
        # inflection point.
        spline2 = UnivariateSpline(temps, derivs, s=3)

        # The root of the 2nd-derivative spline gives the inflection point.
        return spline2.roots()[0], spline

    def alteration_index(self):
        """Return alteration index."""
        return self.heating[1][0] - self.cooling[1][0]

    def correct_for_volume(self, data):
        scale = self.nom_vol / self.real_vol
        return [scale * datum for datum in data]

    @staticmethod
    def shunt_up(values):
        """Ensure that a list of scalars is non-negative.

        If min(values)<0, return values - min(values),
        otherwise return values.

        Args:
          values (list): magnetic suscepetibility values

        Returns:
          values, incremented by a constant 

        """
        if len(values) == 0:
            return values
        minimum = min(values)
        if minimum < 0:
            values = [v - minimum for v in values]
        return values


class MeasurementSet:
    """The results of a series of heating-cooling cycles on a single sample."""

    @staticmethod
    def shunt(heat_cool, offset):
        heat, cool = heat_cool
        heat_s = (heat[0], [m + offset for m in heat[1]])
        cool_s = (cool[0], [m + offset for m in cool[1]])
        return heat_s, cool_s

    def make_zero_at_700(self):
        """Correct values for a zero susceptibility at/near 700 degrees"""
        print(self.name, self.cycles.keys(), self.cycles[700][0][1][:5])
        offset = -min(self.cycles[700][0][1][-5:])
        new_data = {}
        for temp in self.cycles.keys():
            new_data[temp] = MeasurementSet.shunt(self.cycles[temp], offset)
        self.cycles = new_data

    @staticmethod
    def filename_to_temp(filename):
        """Convert a filename to a temperature"""
        leafname = os.path.basename(filename)
        m = re.search(r'^(\d+)[AB]?\.CUR$', leafname)
        if m is None:
            return None
        return int(m.group(1))

    def set_oom(self, new_oom):
        scale = 10. ** (self.oom - new_oom)
        new_data = {}
        for (temp, (heating, cooling)) in self.cycles.items():
            heating2 = (heating[0], [ms*scale for ms in heating[1]])
            cooling2 = (cooling[0], [ms*scale for ms in cooling[1]])
            new_data[temp] = (heating2, cooling2)
        self.cycles = new_data
        self.oom = new_oom

    def read_files(self, sample_dir):
        cur_files = glob.glob(os.path.join(sample_dir, "*.CUR"))
        for filename in cur_files:
            temperature = MeasurementSet.filename_to_temp(filename)
            if temperature is None:
                continue
            self.cycles[temperature] = MeasurementCycle(
                self.furnace, filename, self.real_vol, self.nom_vol)
        # self.make_zero_at_700()

    def __init__(self, furnace, sample_dir, real_vol=0.25, nom_vol=10.):
        self.oom = -6.  # order of magnitude
        self.name = os.path.basename(sample_dir)
        self.furnace = furnace
        self.cycles = {}
        self.nom_vol = nom_vol
        self.real_vol = real_vol
        if sample_dir is not None:
            self.read_files(sample_dir)

    @staticmethod
    def alterations(cycles):
        return [cycle.alteration() for cycle in cycles]
