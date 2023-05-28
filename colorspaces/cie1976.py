"""
Defines L*a*b and L*u*v color space related functions and a color difference formula based on L*a*b
"""

from numpy import asarray, cbrt, divide, einsum, linalg, moveaxis
from numpy.linalg import norm, solve


def transform(a, b):
    return einsum("...ij,...j", a, b, optimize=True)


def solve(a, b):
    return moveaxis(linalg.solve(a, moveaxis(b, -1, 0).reshape(3, -1)), 0, -1).reshape(
        b.shape
    )


D65_XYZ = asarray([0.95047155, 1.0, 1.08882966])

CIE_E = 216 / 24389  # http://www.brucelindbloom.com/LContinuity.html
CIE_K = 24389 / 27  # http://www.brucelindbloom.com/LContinuity.html


def gamma(t1):
    t2 = asarray(cbrt(t1, where=t1 > CIE_E))
    t2[t1 <= CIE_E] = (t1[t1 <= CIE_E] * CIE_K + 16) / 116
    return t2


def gamma_inv(t2):
    t1 = asarray(t2**3)
    t1[t1 <= CIE_E] = (t2[t1 <= CIE_E] * 116 - 16) / CIE_K
    return t1


def _u(X, Y, Z):
    return divide(4 * X, X + 15 * Y + 3 * Z)


def _v(X, Y, Z):
    return divide(9 * Y, X + 15 * Y + 3 * Z)


def Luv_from_XYZ(XYZ, W=D65_XYZ):
    X, Y, Z = moveaxis(XYZ, -1, 0)

    L = 116 * gamma(Y / W[1]) - 16
    u = 13 * L * (_u(X, Y, Z) - _u(*W))
    v = 13 * L * (_v(X, Y, Z) - _v(*W))

    return moveaxis([L, u, v], 0, -1)


def XYZ_from_Luv(Luv, W=D65_XYZ):
    L, u, v = moveaxis(Luv, -1, 0)
    Y = gamma_inv((L + 16) / 116)

    a = (divide(52 * L, u + 13 * L * _u(*W)) - 1) / 3
    b = -5 * Y
    c = -1 / 3
    d = Y * (divide(39 * L, v + 13 * L * _v(*W)) - 5)

    X = divide(d - b, a - c)
    Z = X * a + b
    return moveaxis([X, Y, Z], 0, -1)


_Lab_from_gamma_XYZ = asarray(
    [
        [0, 116, 0],
        [500, -500, 0],
        [0, 200, -200],
    ]
)


def Lab_from_XYZ(XYZ, W=D65_XYZ):
    "http://www.brucelindbloom.com/Eqn_XYZ_to_Lab.html"
    XYZ = asarray(XYZ)
    return transform(_Lab_from_gamma_XYZ, gamma(XYZ / W)) - [16, 0, 0]


def XYZ_from_Lab(Lab, W=D65_XYZ):
    "http://www.brucelindbloom.com/Eqn_XYZ_to_Lab.html"
    Lab = asarray(Lab)
    # solve(_Lab_from_gamma_XYZ, moveaxis(
    #     Lab, -1, 0).reshape(3, -1)).reshape(Lab.shape)
    return gamma_inv(solve(_Lab_from_gamma_XYZ, Lab + [16, 0, 0])) * W


def delta_e_cie1976_from_Lab(Lab1, Lab2):
    "Delta E (CIE 1976) color difference formula for Lab color space"
    return norm(Lab1 - Lab2, axis=-1)
