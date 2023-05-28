"""
Defines CIE Lab color space related functions and delta_e_cie1976 color difference formula
"""

from numpy import asarray, cbrt, einsum, linalg, moveaxis, ndarray
from numpy.linalg import norm
from numpy.typing import ArrayLike


def transform(a, b):
    return einsum("...ij,...j", a, b, optimize=True)


def solve(a, b):
    return moveaxis(linalg.solve(a, moveaxis(b, -1, 0).reshape(3, -1)), 0, -1).reshape(
        b.shape
    )


D65 = {
    "XYZ": asarray([0.95047155, 1.0, 1.08882966]),
}

CIE_E = 216 / 24389  # http://www.brucelindbloom.com/LContinuity.html
CIE_K = 24389 / 27  # http://www.brucelindbloom.com/LContinuity.html


def gamma(t1: ArrayLike) -> ndarray:
    t2 = asarray(cbrt(t1, where=t1 > CIE_E))
    t2[t1 <= CIE_E] = (t1[t1 <= CIE_E] * CIE_K + 16) / 116
    return t2


def gamma_inv(t2: ArrayLike) -> ndarray:
    t1 = asarray(t2**3)
    t1[t1 <= CIE_E] = (t2[t1 <= CIE_E] * 116 - 16) / CIE_K
    return t1


_Lab_from_xyz = asarray(
    [
        [0, 116, 0],
        [500, -500, 0],
        [0, 200, -200],
    ]
)


def Lab_from_XYZ(XYZ: ArrayLike, W=D65["XYZ"]) -> ndarray:
    "http://www.brucelindbloom.com/Eqn_XYZ_to_Lab.html"
    XYZ = asarray(XYZ)
    return transform(_Lab_from_xyz, gamma(XYZ / W)) - [16, 0, 0]


def XYZ_from_Lab(Lab: ArrayLike, W=D65["XYZ"]) -> ndarray:
    "http://www.brucelindbloom.com/Eqn_XYZ_to_Lab.html"
    Lab = asarray(Lab)
    return gamma_inv(solve(_Lab_from_xyz, Lab + [16, 0, 0])) * W


def delta_e_cie1976_from_Lab(Lab1, Lab2):
    "Delta E (CIE 1976) color difference formula for Lab color space"
    return norm(Lab1 - Lab2, axis=-1)
