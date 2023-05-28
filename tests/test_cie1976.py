from colorspaces import cie1976
from conftest import _is_inverse, _test_delta


def test_gamma():
    _is_inverse(cie1976.gamma, cie1976.gamma_inv)


def test_Luv_XYZ():
    _is_inverse(cie1976.Luv_from_XYZ, cie1976.XYZ_from_Luv)


def test_Lab_XYZ():
    _is_inverse(cie1976.Lab_from_XYZ, cie1976.XYZ_from_Lab)


def test_delta_e_cie1976_from_Lab():
    _test_delta(cie1976.delta_e_cie1976_from_Lab, 0.37416573867739417)
