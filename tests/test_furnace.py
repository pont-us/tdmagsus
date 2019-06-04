#!/usr/bin/python3

from tdmagsus import Furnace
from numpy import array, array_equal

import unittest


class TestFurnace(unittest.TestCase):

    def test_extend_data(self):
        temps = array([30, 40, 50, 60])
        mss = array([21, 22, 23, 24])
        temps_result, mss_result = Furnace.extend_data((temps, mss))
        self.assertTrue(array_equal(array([10, 20, 30, 40, 50, 60, 70, 80]),
                                    temps_result))
        self.assertTrue(array_equal(array([21, 21, 21, 22, 23, 24, 24, 24]),
                                    mss_result))

    def test_correct_with_spline(self):
        def spline(temperature: float) -> float:
            spline_dict = {30.: 5, 40.: 4, 50.: 3, 60.: 2}
            return spline_dict[temperature]

        temps = [30., 40., 50., 60.]
        mss = [21, 22, 23, 24]

        temps_result, mss_result = \
            Furnace.correct_with_spline(temps, mss, spline)

        self.assertListEqual(temps, temps_result)
        self.assertTrue(array_equal(array([16, 18, 20, 22]), mss_result))
