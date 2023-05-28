from colorspaces import cie2000
from conftest import _test_delta


def test_delta_e_cie2000_from_Lab():
    _test_delta(cie2000.delta_e_cie2000_from_Lab, 0.4164154127016923)


# def test_delta_e_cmc_from_LCh():
#     _test_delta(cie2000.delta_e_cmc_from_LCh, 17.063221085266047)


# def test_delta_e_cmc_from_Lab():
#     _test_delta(cie2000.delta_e_cmc_from_Lab, 15.716265438196919)
