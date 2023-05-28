"""
Defines sRGB, XYZ and xyY color space related functions
"""

from numpy import (
    array,
    asarray,
    divide,
    einsum,
    empty_like,
    linalg,
    moveaxis,
    ndarray,
    power,
    sum,
)
from numpy.linalg import solve
from numpy.typing import ArrayLike

# D65 White point in xyY color space
D65 = {"xyY": array([0.312727, 0.329023, 1])}

# sRGB triangle in xyY color space
sRGB = {
    "xyY": array(
        [
            [0.64, 0.33, 1],  # red
            [0.30, 0.60, 1],  # green
            [0.15, 0.06, 1],  # blue
        ]
    )
    # array([
    #     [.64, .33, .21267395],  # red
    #     [.30, .60, .71515107],  # green
    #     [.15, .06, .07217497],  # blue
    # ])
}


def transform(a: ArrayLike, b: ArrayLike) -> ndarray:
    return einsum("...ij,...j", a, b, optimize=True)


def solve(a, b):
    return moveaxis(linalg.solve(a, moveaxis(b, -1, 0).reshape(3, -1)), 0, -1).reshape(
        b.shape
    )


def XYZ_from_xyY(xyY: ArrayLike) -> ndarray:
    "http://www.brucelindbloom.com/Eqn_xyY_to_XYZ.html"

    x, y, Y = moveaxis(xyY, -1, 0)

    X = asarray(divide(x * Y, y, where=y != 0))
    Z = asarray(divide((1 - x - y) * Y, y, where=y != 0))

    X[y == 0] = 0
    Z[y == 0] = 0

    return moveaxis([X, Y, Z], 0, -1)


def xyY_from_XYZ(XYZ: ArrayLike) -> ndarray:
    "http://www.brucelindbloom.com/Eqn_XYZ_to_xyY.html"

    S = sum(XYZ, axis=-1)
    # xyz = divide(XYZ, S, where=S != 0)
    # xyz[S == 0] = 0

    X, Y, Z = moveaxis(XYZ, -1, 0)

    x = asarray(divide(X, S, where=S != 0))
    y = asarray(divide(Y, S, where=S != 0))

    x[S == 0] = D65["xyY"][0]  # .312727
    y[S == 0] = D65["xyY"][1]  # .329023

    return moveaxis([x, y, Y], 0, -1)


def _calc_matrix_XYZ_from_RGB(
    RGB=XYZ_from_xyY(sRGB["xyY"]).T, W=XYZ_from_xyY(D65["xyY"])
):
    "http://www.brucelindbloom.com/Eqn_RGB_XYZ_Matrix.html"

    S = solve(RGB, W)
    return S * RGB


_XYZ_from_RGB = _calc_matrix_XYZ_from_RGB()

# _XYZ_from_RGB = array([
#     [0.41245858, 0.35757554, 0.18043743],
#     [0.21267395, 0.71515107, 0.07217497],
#     [0.01933400, 0.11919185, 0.95030382],
# ])

# M = array([
#     [0.41245643, 0.35757608, 0.18043748],
#     [0.21267285, 0.71515217, 0.072175],
#     [0.0193339, 0.11919203, 0.95030407],
# ])

# M_=array([
#     [3.2404542, -1.5371385, -0.4985314],
#     [-0.9692660, 1.8760108, 0.0415560],
#     [0.0556434, -0.2040259, 1.0572252],
# ]),


def XYZ_from_RGB(RGB: ArrayLike) -> ndarray:
    "convert linear rgb values to XYZ values"
    return transform(_XYZ_from_RGB, RGB)


def RGB_from_XYZ(XYZ: ArrayLike) -> ndarray:
    return solve(_XYZ_from_RGB, XYZ)


def RGB_from_sRGB(rgb: ArrayLike, gamma="sRGB") -> ndarray:
    "Gamma correction"

    rgb = asarray(rgb)
    if gamma == "sRGB":
        RGB = power((rgb + 0.055) / 1.055, 2.4, where=rgb > 0.04045)
        RGB[rgb <= 0.04045] = rgb[rgb <= 0.04045] / 12.92
        return RGB
    else:
        return power(rgb, gamma)


def sRGB_from_RGB(RGB: ArrayLike, gamma="sRGB") -> ndarray:
    "Gamma correction"

    if gamma == "sRGB":
        rgb = power(RGB, 1 / 2.4, where=RGB > 0.0031308) * 1.055 - 0.055
        rgb[RGB <= 0.0031308] = RGB[RGB <= 0.0031308] * 12.92
        return rgb
    else:
        return power(RGB, 1 / gamma)


def sRGB_from_XYZ(XYZ: ArrayLike) -> ndarray:
    return sRGB_from_RGB(RGB_from_XYZ(XYZ))


def XYZ_from_sRGB(rgb: ArrayLike) -> ndarray:
    return XYZ_from_RGB(RGB_from_sRGB(rgb))


def xyY_from_T(T: ArrayLike) -> ndarray:
    x = empty_like(T, dtype=float)
    x[T <= 7000] = -4.607e9 / T**3 + 2.9678e6 / T**2 + 0.09911e3 / T + 0.244063
    x[T > 7000] = -2.0064e9 / T**3 + 1.9018e6 / T**2 + 0.24748e3 / T + 0.237040

    y = -3 * x**2 + 2.87 * x - 0.275

    return moveaxis([x, y, 1], 0, -1)
