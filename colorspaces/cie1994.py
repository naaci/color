"""
Defines CIE 1994 CIE LCHab color space and color difference formula
"""

from numpy import arctan2, asarray, cos, divide, moveaxis, pi, sin, sqrt


def LCh_from_Lab(Lab):
    L, a, b = moveaxis(Lab, -1, 0)
    C = sqrt(a**2 + b**2)
    h = asarray(arctan2(b, a)) % (2 * pi)

    return moveaxis([L, C, h], 0, -1)


def Lab_from_LCh(LCh):
    L, C, h = moveaxis(LCh, -1, 0)
    a = C * cos(h)
    b = C * sin(h)

    return moveaxis([L, a, b], 0, -1)


def LCh_from_Luv(Luv):
    L, u, v = moveaxis(Luv, -1, 0)
    C = sqrt(u**2 + v**2)
    h = asarray(arctan2(v, u)) % (2 * pi)

    return moveaxis([L, C, h], 0, -1)


def Luv_from_LCh(LCh):
    L, C, h = moveaxis(LCh, -1, 0)
    u = C * cos(h)
    v = C * sin(h)

    return moveaxis([L, u, v], 0, -1)


def delta_e_cie1994_from_LCh(LCh1, LCh2, textile=False):
    "Delta E (CIE 1994) color difference formula for LChab color space"

    L1, C1, h1 = moveaxis(LCh1, -1, 0)
    L2, C2, h2 = moveaxis(LCh2, -1, 0)

    ΔL = L1 - L2
    ΔC = C1 - C2
    Δh = h1 - h2

    if textile:
        K1 = 0.048
        K2 = 0.014
        K_L = 2
    else:
        K1 = 0.045
        K2 = 0.015
        K_L = 1

    K_C = 1
    K_h = 1

    S_L = 1
    S_C = 1 + K1 * C1
    S_h = 1 + K2 * C2

    return sqrt(
        divide(ΔL, K_L * S_L) ** 2
        + divide(ΔC, K_C * S_C) ** 2
        + abs(divide(Δh, K_h * S_h))
    )


def delta_e_cie1994_from_Lab(Lab1, Lab2):
    "Delta E (CIE 1994) color difference formula for L*a*b color space"

    return delta_e_cie1994_from_LCh(
        LCh_from_Lab(Lab1),
        LCh_from_Lab(Lab2),
    )
