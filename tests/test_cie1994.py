from cie1994 import *
from conftest import _is_inverse, _test_delta


def test_LCh_Lab_LCh():
    _is_inverse(LCh_from_Lab, Lab_from_LCh)


def test_Lab_LCh_Lab():
    _is_inverse(Lab_from_LCh, LCh_from_Lab)


def test_delta_e_cie1994_from_Lab():
    _test_delta(delta_e_cie1994_from_Lab, 0.368621751467173)
