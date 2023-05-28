from colorspaces import cie1931
from conftest import _is_inverse


def test_XYZ_xyY():
    _is_inverse(cie1931.XYZ_from_xyY, cie1931.xyY_from_XYZ)


def test_XYZ_RGB():
    _is_inverse(cie1931.XYZ_from_RGB, cie1931.RGB_from_XYZ)


def test_sRGB_RGB():
    _is_inverse(cie1931.RGB_from_sRGB, cie1931.sRGB_from_RGB)
