#!/usr/bin/python3

from tdmagsus import read_cur_file
import os.path
import hashlib

import unittest


class TestFunctions(unittest.TestCase):
    """Tests for functions outside classes"""

    def test_read_cur_file(self):
        filename = os.path.join(os.path.dirname(__file__), "testdata",
                                "TUBE1.CUR")
        (heating_t, heating_ms), (cooling_t, cooling_ms) = \
            read_cur_file(filename)

        def check_hash(expected, array):
            # We use hashlib here, because since Python 3.3 object.__hash__
            # incorporates a random salt, producing different hash values
            # between program runs for the same input. See
            # https://docs.python.org/3/reference/datamodel.html#object.__hash__
            self.assertEqual(expected,
                             hashlib.sha224(array.tobytes()).hexdigest())

        # Characterization tests with precalculated hashes for the arrays --
        # more concise and convenient than using full arrays for the expected
        # values.
        check_hash("1ac194f68a365c9b7524d7f9146487a685a7fd4a999764ddec174147",
                   heating_t)
        check_hash("ac19cce10c9380be6ccc761746488a5d189ec2d41cfeb6df6d755abf",
                   heating_ms)
        check_hash("07d78fda1292e48c63f58ecb46a7c7e98c01189090ee7fc27bf995ce",
                   cooling_t)
        check_hash("999d9229f8ecaae3d56b5a3b9d490bb8328ea50f12260928de5ab6fd",
                   cooling_ms)
