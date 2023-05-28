from hunterlab import *
from conftest import _is_inverse


def test_XYZ_Lab():
    _is_inverse(XYZ_from_Lab, Lab_from_XYZ)
