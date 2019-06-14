#!/usr/bin/python3

from tdmagsus import MeasurementSet

import unittest


class TestMeasurementSet(unittest.TestCase):

    def test_shunt(self):
        (heat_t, heat_ms), (cool_t, cool_ms) = \
            MeasurementSet.shunt(
                (([1, 2, 3, 4], [10, 20, 30, 40]),
                 ([4, 3, 2, 1], [45, 35, 25, 15])),
                2)
        self.assertListEqual([1, 2, 3, 4], heat_t)
        self.assertListEqual([4, 3, 2, 1], cool_t)
        self.assertListEqual([12, 22, 32, 42], heat_ms)
        self.assertListEqual([47, 37, 27, 17], cool_ms)

    def test_filename_to_temp(self):
        self.assertIsNone(MeasurementSet.filename_to_temp("unparseable"))
        self.assertEqual(700, MeasurementSet.filename_to_temp("700.CUR"))
        self.assertEqual(50, MeasurementSet.filename_to_temp("50A.CUR"))
        self.assertEqual(5, MeasurementSet.filename_to_temp("5B.CUR"))
