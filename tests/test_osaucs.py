from colorspaces import osaucs
from numpy import asarray, isclose


def test_Lgj_from_XYZ():
    c = asarray([0.1, 0.2, 0.3])
    assert isclose(
        osaucs.Lgj_from_XYZ(c), [-10.99395228, 5.68710901, -1.29335382]
    ).all()
