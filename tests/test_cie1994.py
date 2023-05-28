from colorspaces import cie1994
from conftest import _is_inverse, _test_delta


def test_LCh_Lab_LCh():
    _is_inverse(cie1994.LCh_from_Lab, cie1994.Lab_from_LCh)


def test_Lab_LCh_Lab():
    _is_inverse(cie1994.Lab_from_LCh, cie1994.LCh_from_Lab)


def test_delta_e_cie1994_from_Lab():
    _test_delta(cie1994.delta_e_cie1994_from_Lab, 0.368621751467173)
