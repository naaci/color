from cie1931 import *
from conftest import _is_inverse


def test_XYZ_xyY():
    _is_inverse(XYZ_from_xyY, xyY_from_XYZ)


def test_XYZ_RGB():
    _is_inverse(XYZ_from_RGB, RGB_from_XYZ)


def test_sRGB_RGB():
    _is_inverse(RGB_from_sRGB, sRGB_from_RGB)
