"""
A perceptual color space for image processing

https://bottosson.github.io/posts/oklab/
"""

from numpy import array, cbrt, einsum, linalg, moveaxis, ndarray
from numpy.linalg import solve
from numpy.typing import ArrayLike


def transform(a, b):
    return einsum("...ij,...j", a, b, optimize="greedy")


def solve(a, b):
    return moveaxis(linalg.solve(a, moveaxis(b, -1, 0).reshape(3, -1)), 0, -1).reshape(
        b.shape
    )


_LMS_from_XYZ = array(
    [
        [0.8189330101, 0.3618667424, -0.1288597137],
        [0.0329845436, 0.9293118715, 0.0361456387],
        [0.0482003018, 0.2643662691, 0.6338517070],
    ]
)

_Lab_from_lms = array(
    [
        [0.2104542553, 0.7936177850, -0.0040720468],
        [1.9779984951, -2.4285922050, 0.4505937099],
        [0.0259040371, 0.7827717662, -0.8086757660],
    ]
)


def LMS_from_XYZ(XYZ: ArrayLike) -> ndarray:
    return transform(_LMS_from_XYZ, XYZ)


def XYZ_from_LMS(SML: ArrayLike) -> ndarray:
    return solve(_LMS_from_XYZ, SML)


def Lab_from_LMS(SML: ArrayLike) -> ndarray:
    return transform(_Lab_from_lms, cbrt(SML))


def LMS_from_Lab(LAB: ArrayLike) -> ndarray:
    return solve(_Lab_from_lms, LAB) ** 3


def Lab_from_XYZ(XYZ: ArrayLike) -> ndarray:
    return Lab_from_LMS(LMS_from_XYZ(XYZ))


def XYZ_from_Lab(LAB: ArrayLike) -> ndarray:
    return XYZ_from_LMS(LMS_from_Lab(LAB))
