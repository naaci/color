from oklab import *
from conftest import _is_inverse, _test_delta


def test_LMS_Lab():
    _is_inverse(LMS_from_Lab, Lab_from_LMS)


def test_XYZ_LMS():
    _is_inverse(XYZ_from_LMS, LMS_from_XYZ)
