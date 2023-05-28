from colorspaces import cielab
from conftest import _is_inverse, _test_delta


def test_gamma():
    _is_inverse(cielab.gamma, cielab.gamma_inv)


def test_Lab_XYZ():
    _is_inverse(cielab.Lab_from_XYZ, cielab.XYZ_from_Lab)


def test_delta_e_cie1976_from_Lab():
    _test_delta(cielab.delta_e_cie1976_from_Lab, 0.37416573867739417)
