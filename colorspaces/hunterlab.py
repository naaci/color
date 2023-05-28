"""
Hunter Lab
"""

from numpy import array, asarray, divide, moveaxis, ndarray, sqrt
from numpy.typing import ArrayLike

D65 = {
    "K": array([100, 172.3, 67.2]),
    "XYZ": asarray([0.95047155, 1.0, 1.08882966]),
}


def Lab_from_XYZ(XYZ: ArrayLike, W=D65) -> ndarray:
    XYZ = asarray(XYZ)
    x, y, z = moveaxis(XYZ / W["XYZ"], -1, 0)

    L = sqrt(y)
    a = divide(x - y, L)
    b = divide(y - z, L)

    return moveaxis([L, a, b], 0, -1) * W["K"]


def XYZ_from_Lab(Lab: ArrayLike, W=D65) -> ndarray:
    Lab = asarray(Lab)
    L, a, b = moveaxis(Lab / W["K"], -1, 0)

    y = L**2
    x = y + L * a
    z = y - L * b

    return moveaxis([x, y, z], 0, -1) * W["XYZ"]
