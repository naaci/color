from colorspaces import oklab
from conftest import _is_inverse


def test_LMS_Lab():
    _is_inverse(oklab.LMS_from_Lab, oklab.Lab_from_LMS)


def test_XYZ_LMS():
    _is_inverse(oklab.XYZ_from_LMS, oklab.LMS_from_XYZ)
