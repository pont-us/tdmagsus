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

"""Tools for working with temperature-dependent magnetic susceptibility data.

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
from typing import List, Optional, Tuple, Callable

import numpy
from numpy import array, ndarray
from scipy.interpolate import UnivariateSpline

line_pattern = re.compile(r"^ +\d")
field_separator = re.compile(" +")


def read_cur_file(filename: str) ->\
        Tuple[Tuple[ndarray, ndarray], Tuple[ndarray, ndarray]]:
    """Read a .CUR magnetic susceptibility file.

    :param filename: name of file to read
    :return: ((heating_temps, heating_mag_sus_values),
       (cooling_temps, cooling_mag_sus_values))
      This is a tuple of two tuples, each containing two ndarrays.
    """

    temperatures = []
    mag_suss = []
    with open(filename, "r") as infile:
        for line in infile:
            if line_pattern.match(line.rstrip()):
                temperature, mag_sus = \
                    map(float, field_separator.split(line.lstrip())[0:2])
                temperatures.append(temperature)
                mag_suss.append(mag_sus)

    # Looking for the first temperature decrease is an unreliable way to find
    # the start of the cooling curve: it's not guaranteed that the heating
    # curve will be monotonically increasing, especially around 100°C where
    # evaporation from a moist sample can cause a brief drop in temperature.
    # Instead we split at the maximum of the entire temperature series.
    split_at = temperatures.index(max(temperatures))
    # The maximum temperature is duplicated: it's the last step in the
    # heating array, and also the first step in the cooling array.
    heating = (temperatures[:(split_at + 1)], mag_suss[:(split_at + 1)])
    cooling = (temperatures[split_at:], mag_suss[split_at:])

    cooling[0].reverse()
    cooling[1].reverse()
    heating = (array(heating[0][1:]), array(heating[1][1:]))
    cooling = (array(cooling[0][1:]), array(cooling[1][1:]))
    return heating, cooling


class Furnace:
    """The temperature-susceptibility behaviour of an empty furnace.

    This class represents the susceptibility of an empty furnace during a
    heating-cooling cycle. It is used to correct the measured susceptibility
    of a sample by substracting the furnace's susceptibility.
    """

    @staticmethod
    def extend_data(temps_mss: Tuple[ndarray, ndarray]) \
            -> Tuple[ndarray, ndarray]:
        """Pad data at ends to improve spline fit.

        Given ([T1, T2, T3, ... , Tn-1, Tn], [M1, M2, M3, ... , Mn-1, Mn]),
        produces ([T1-20, T1-10, T1, T2, T3, ... , Tn-1, Tn, Tn+10, Tn+20],
                  [M1, M1, M1, M2, M3, ... , Mn-1, Mn, Mn, Mn]).

        :param temps_mss: tuple of temperatures and mag sus values
        :return: same values, padded by two extra data points at each end
        """
        temps, mss = temps_mss
        tlist = list(temps.tolist())
        mlist = list(mss.tolist())
        tlist = [(tlist[0] - 20)] + [(tlist[0] - 10)] + \
            tlist + [tlist[-1] + 10] + [tlist[-1] + 20]
        mlist = [mlist[0]] + [mlist[0]] + \
            mlist + [mlist[-1]] + [mlist[-1]]
        return array(tlist), array(mlist)

    def __init__(self, filename: str, smoothing: float = 100) -> None:
        """Initialize Furnace object from a CUR file.

        :param filename: file path from which to read furnace data
        :param smoothing: smoothing factor for spline curve
        """
        heat, cool = read_cur_file(filename)
        self.heat_data = Furnace.extend_data(heat)
        self.cool_data = Furnace.extend_data(cool)
        self.heat_spline = UnivariateSpline(*self.heat_data, s=smoothing)
        self.cool_spline = UnivariateSpline(*self.cool_data, s=smoothing)

    def get_spline_data(self) -> \
            Tuple[Tuple[ndarray, ndarray], Tuple[ndarray, ndarray],
                  Tuple[ndarray, ndarray], Tuple[ndarray, ndarray]]:
        """Return furnace temperature/M.S. data and spline approximations.
        
        This method is mainly intended for checking that the splines are
        doing a good job of smoothing the data.

        :return: (heating_data, heating_spline, cooling_data, cooling_spline).
           Each element of this tuple is itself a 2-tuple containing an
           ndarray of temperatures and an ndarray of associated M.S. values.
        """

        spline_x = numpy.arange(20, 701)
        spline_y_heating = self.heat_spline(spline_x)
        spline_y_cooling = self.cool_spline(spline_x)
        return self.heat_data, (spline_x, spline_y_heating),\
            self.cool_data, (spline_x, spline_y_cooling)

    @staticmethod
    def correct_with_spline(temps: List[float], mss: List[float],
                            spline: Callable[[float], float])\
            -> Tuple[List[float], ndarray]:
        """Correct susceptibility measurements using a supplied spline.

        :param temps: an array_like of temperatures
        :param mss: an array_like of magnetic susceptibility measurement
                    taken at the specified temperatures
        :param spline: a function mapping a temperature to a furnace
                       susceptibility
        :return: a 2-tuple containing temperatures and a corresponding series of
                 corrected susceptibilities
        """
        mss_corrected = numpy.zeros_like(mss)
        for i in range(0, len(temps)):
            mss_corrected[i] = mss[i] - spline(temps[i])
        return temps, mss_corrected

    def correct(self, heating: Tuple, cooling: Tuple) -> Tuple[Tuple, Tuple]:
        """Correct susceptibility measurements using these furnace measurements

        :param heating: a 2-tuple (temperatures, susceptibilities) of
                        array_likes for the heating measurements
        :param cooling: a 2-tuple (temperatures, susceptibilities) of
                        array_likes for the cooling measurements
        :return: a 2-tuple of 2-tuples giving the corrected measurements:
                 ((heating_temperatures, heating_susceptibilities),
                  (cooling_temperatures, cooling_susceptibilities))
        """
        return (Furnace.correct_with_spline(heating[0], heating[1],
                                            self.heat_spline),
                Furnace.correct_with_spline(cooling[0], cooling[1],
                                            self.cool_spline))


class MeasurementCycle:
    """The results of a single heating-cooling run."""

    def __init__(self, furnace: Furnace, filename: str,
                 real_vol: float = 0.25, nom_vol: float = 10.0) -> None:
        """Create a measurement cycle object from data in a specified file.

        :param furnace: empty furnace correction
        :param filename:
        :param real_vol:
        :param nom_vol:
        """

        self.furnace = furnace
        self.real_vol = real_vol
        self.nom_vol = nom_vol
        (heating, cooling) = read_cur_file(filename)
        if self.furnace is not None:
            heating, cooling = self.furnace.correct(heating, cooling)
        heating = (heating[0], self.correct_for_volume(heating[1]))
        cooling = (cooling[0], self.correct_for_volume(cooling[1]))
        self.data = (heating, cooling)

    def write_csv(self, filename: str) -> None:
        """Write furnace-corrected data to a CSV file.

        :param filename: name of file to write to.
        """

        cooling = (list(reversed(self.data[1][0])),
                   list(reversed(self.data[1][1])))
        with open(filename, "w") as fh:
            for direction in self.data[0], cooling:
                pairs = zip(direction[0], direction[1])
                for pair in pairs:
                    fh.write("%.2f,%.2f\n" % (pair[0], pair[1]))

    @staticmethod
    def chop_data(temps_mss: Tuple[ndarray, ndarray], min_temp: float,
                  max_temp: float) -> Tuple[ndarray, ndarray]:
        """Truncate data to a given temperature range.

        Any temperature falling outside the specified range will be excluded
        from the returned temperature array, and its corresponding
        susceptibility value will be excluded from the returned susceptibility
        array.

        :param temps_mss: a 2-tuple of lists, (temperatures, susceptibilites)
        :param min_temp: minimum temperature for truncation
        :param max_temp: maximum remperature for truncation
        :return: a 2-tuple of ndarrays, (temperatures, susceptibilities), where
                 all temperatures are in the range (min_temp, max_temp)
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
    def linear_fit(xs: ndarray, ys: ndarray) -> Tuple[numpy.poly1d, float]:
        """Make a least-squares linear fit to a 2-dimensional data set.

        This method also calculates and returns the r-squared value associated
        with the fit.

        :param xs: x co-ordinates
        :param ys: y co-ordinates
        :return: a tuple of (numpy.poly1d, r_squared)
        """
        fit = numpy.polyfit(xs, ys, 1)
        poly = numpy.poly1d(fit.tolist())
        model_ys = poly(xs)
        mean_y = numpy.mean(ys)
        sserr = numpy.sum((ys - model_ys) ** 2)
        sstot = numpy.sum((ys - mean_y) ** 2)
        rsquared = 1 - sserr / sstot
        # rsquared is already a float; the "conversion" is to help type checkers
        return poly, float(rsquared)

    def curie_paramag(self, min_temp: float, max_temp: float)\
            -> Tuple[float, float, numpy.poly1d]:
        """Estimate Curie temperature by linear fit to inverse susceptibility.

        :param min_temp: minimum of temperature range for fit
        :param max_temp: maximum of temperature range for fit
        :return: (curie, rsquared, poly) where
          curie is estimated Curie temperature;
          rsquared is R² value for fit;
          poly is polynomial object representing line of best fit.
        """

        temps, mss = \
            MeasurementCycle.chop_data(self.data[0], min_temp, max_temp)
        poly, rsquared = MeasurementCycle.linear_fit(temps, 1. / mss)
        curie = poly.r[0]  # x axis intercept
        return curie, rsquared, poly

    def curie_inflection(self, min_temp: float, max_temp: float)\
            -> Tuple[float, UnivariateSpline]:
        """Estimate Curie temperature by inflection point.
        
        Estimate Curie point by determining the inflection point of
        the curve segment starting at the Hopkinson peak. The curve
        segment must be specified.

        :param min_temp: start of curve segment
        :param max_temp: end of curve segment
        :return: (temp, spline) where
          temp is the estimated Curie temperature;
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

    def alteration_index(self) -> float:
        """Return the alteration index for this cycle.

        :return: the alteration idex for this cycle
        """
        return self.data[0][1][0] - self.data[1][1][0]

    def correct_for_volume(self, data: ndarray) -> ndarray:
        """Correct supplied data for volume.

        The volume correction factor is (nominal_volume / real_volume). The
        values are set by the constructor.

        :param data: a list of numerical data to correct
        :return: a list containing the supplied data scaled by the
                 volume correction factor
        """
        scale = self.nom_vol / self.real_vol
        return array([scale * datum for datum in data])

    @staticmethod
    def shunt_up(values: List[float]) -> List[float]:
        """Ensure that a list of scalars is non-negative.
        
        If min(values) < 0, return values - min(values),
        otherwise return values.

        :param values: magnetic suscepetibility values
        :type values: list
        :return: values, incremented by a constant
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
    def shunt(heat_cool: Tuple[Tuple[List[float], List[float]],
                               Tuple[List[float], List[float]]], offset: float)\
            -> Tuple[Tuple[List[float], List[float]],
                     Tuple[List[float], List[float]]]:
        """Offset all magnetic susceptibility values by a supplied value.

        :param heat_cool: ((heating_temps, heating_susceptibilities),
                           (cooling_temps, cooling_susceptibilities))
        :param offset: amount to add to or subtract from each susceptibility
        :return: a tuple of tuples like heat_cool, but with susceptibilies
                 offset
        """
        heat, cool = heat_cool
        heat_s = (heat[0], [m + offset for m in heat[1]])
        cool_s = (cool[0], [m + offset for m in cool[1]])
        return heat_s, cool_s

    def make_zero_at_700(self) -> None:
        """Correct values for a zero susceptibility at/near 700 degrees"""

        print(self.name, self.cycles.keys(), self.cycles[700][0][1][:5])
        offset = -min(self.cycles[700][0][1][-5:])
        new_data = {}
        for temp in self.cycles.keys():
            new_data[temp] = MeasurementSet.shunt(self.cycles[temp], offset)
        self.cycles = new_data

    @staticmethod
    def filename_to_temp(filename: str) -> Optional[int]:
        """Convert a filename to a temperature.

        Extracts a temperature from a numeric filename with a .CUR suffix.
        The temperature may also be suffixed by A or B.
        e.g. "700.CUR" -> 700, "350B.CUR" -> 350

        :param filename: a filename containing a temperature
        :return: the temperature represented by the filename, or None if the
          filename does not represent a temperature
        """
        leafname = os.path.basename(filename)
        m = re.search(r"^(\d+)[AB]?\.CUR$", leafname)
        if m is None:
            return None
        return int(m.group(1))

    def set_oom(self, new_oom: int) -> None:
        """Change this measurement set's order of magnitude.

        :param new_oom: the new order of magnitude
        """
        scale = 10. ** (self.oom - new_oom)
        new_data = {}
        for (temp, (heating, cooling)) in self.cycles.items():
            heating2 = (heating[0], [ms*scale for ms in heating[1]])
            cooling2 = (cooling[0], [ms*scale for ms in cooling[1]])
            new_data[temp] = (heating2, cooling2)
        self.cycles = new_data
        self.oom = new_oom

    def read_files(self, sample_dir: str) -> None:
        """Read a directory of files into this measurement set.

        :param sample_dir: path of directory containing CUR files
        """
        cur_files = glob.glob(os.path.join(sample_dir, "*.CUR"))
        for filename in cur_files:
            temperature = MeasurementSet.filename_to_temp(filename)
            if temperature is None:
                continue
            self.cycles[temperature] = MeasurementCycle(
                self.furnace, filename, self.real_vol, self.nom_vol)

    def __init__(self, furnace: Furnace, sample_dir: str,
                 real_vol: float = 0.25, nom_vol: float = 10.0) -> None:
        """Create a measurement set from data in a specified directory.

        :param furnace: empty furnace data
        :param sample_dir: directory from which to read samples
        :param real_vol: real sample volume (cm³)
        :param nom_vol: nominal sample volume (cm³)
        """
        self.oom = -6.  # order of magnitude
        self.name = os.path.basename(sample_dir)
        self.furnace = furnace
        self.cycles = {}
        self.nom_vol = nom_vol
        self.real_vol = real_vol
        if sample_dir is not None:
            self.read_files(sample_dir)

    @staticmethod
    def alterations(cycles: List[MeasurementCycle]) -> List[float]:
        """Calculate alteration indices for a list of cycles.

        :param cycles: a list of MeasurementCycles
        :return: a corresponding list of alteration indices
        """
        return [cycle.alteration_index() for cycle in cycles]
