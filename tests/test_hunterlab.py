from colorspaces import hunterlab
from conftest import _is_inverse


def test_XYZ_Lab():
    _is_inverse(hunterlab.XYZ_from_Lab, hunterlab.Lab_from_XYZ)
