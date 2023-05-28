"""
Optical Society of America Uniform Color Scales (OSA-UCS)
"""

from . import cie1931
from numpy import array, cbrt, divide, einsum, moveaxis, ndarray, sqrt
from numpy.linalg import norm
from numpy.typing import ArrayLike


def transform(a, b):
    return einsum("...ij,...j", a, b, optimize=True)


_RGB_from_XYZ = array(
    [
        [0.799, 0.4194, -0.1648],
        [-0.4493, 1.3265, 0.0927],
        [-0.1149, 0.3394, 0.717],
    ]
)

_ab_from_rgb = array(
    [
        [-13.7, 17.7, -4],
        [1.7, 8, -9.7],
    ]
)


def Lgj_from_XYZ(XYZ: ArrayLike) -> ndarray:
    x, y, Y = moveaxis(cie1931.xyY_from_XYZ(XYZ), -1, 0)

    K = (
        4.4934 * x**2
        + 4.3034 * y**2
        - 4.276 * x * y
        - 1.3744 * x
        - 2.5643 * y
        + 1.8103
    )
    Y0 = K * Y
    L = 5.9 * (cbrt(Y0) - 2 / 3 + 0.042 * cbrt(Y0 - 30))
    C = divide(L, 5.9 * (cbrt(Y0) - 2 / 3))

    gj = C * transform(_ab_from_rgb, cbrt(transform(_RGB_from_XYZ, XYZ)))

    return moveaxis([sqrt(1 / 2) * (L - 14.3993), *moveaxis(gj, -1, 0)], 0, -1)


def delta(Lgj1, Lgj2):
    return norm(Lgj1 - Lgj2)
