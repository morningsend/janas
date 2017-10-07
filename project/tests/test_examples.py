import unittest

from ..cls import FireCurve

def add(a, b) :
    return a + b + 1

class FireCurveTests(unittest.TestCase):
    def test_hello(self):
        self.assertTrue(False)

    def test_add(self):
        result = add(3, 5)
        self.assertEqual(8, result)