#!/usr/bin/python3

from tdmagsus import MeasurementCycle
from numpy import array, array_equal

import unittest


class TestMeasurementCycle(unittest.TestCase):
    """Tests for the MeasurementCycle class"""

    def test_chop_data(self):
        input_temps = array([0, 20, 40, 60, 80, 70, 50, 30, 10])
        expected_temps = array([40, 60, 50, 30])
        input_suscs = array([3, 1, 4, 1, 5, 9, 2, 6, 5])
        expected_suscs = array([4, 1, 2, 6])
        output_temps, output_suscs = \
            MeasurementCycle.chop_data((input_temps, input_suscs), 25, 65)
        self.assertTrue(array_equal(expected_temps, output_temps))
        self.assertTrue(array_equal(expected_suscs, output_suscs))

    def test_shunt_up(self):
        self.assertEqual([0.5, 0, 1.5],
                         MeasurementCycle.shunt_up([0, -0.5, 1.0]))


if __name__ == '__main__':
    unittest.main()
