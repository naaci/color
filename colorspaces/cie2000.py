"""
Defines CIE 2000 color difference formula
"""

from numpy import abs, arctan2, asarray, cos, divide, exp, moveaxis, pi, sin, sqrt


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


def delta_e_cmc_from_Lab(Lab1, Lab2, l=1, c=1):
    L1, a1, b1 = moveaxis(Lab1, -1, 0)
    L2, a2, b2 = moveaxis(Lab2, -1, 0)

    L1, C1, H1 = moveaxis(LCh_from_Lab(Lab1), -1, 0)
    L2, C2, H2 = moveaxis(LCh_from_Lab(Lab2), -1, 0)

    return delta_e_cmc(L1, L2, C1, C2, a1, a2, b1, b2, l=l, c=c)


def delta_e_cmc_from_LCh(LCh1, LCh2, l=1, c=1):
    L1, C1, H1 = moveaxis(LCh1, -1, 0)
    L2, C2, H2 = moveaxis(LCh2, -1, 0)

    L1, a1, b1 = moveaxis(Lab_from_LCh(LCh1), -1, 0)
    L2, a2, b2 = moveaxis(Lab_from_LCh(LCh2), -1, 0)

    return delta_e_cmc(L1, L2, C1, C2, a1, a2, b1, b2, l=l, c=c)


def delta_e_cmc(L1, L2, C1, C2, a1, a2, b1, b2, l=1, c=1):
    ΔC = C2 - C1

    ΔL = L1 - L2
    Δa = a1 - a2
    Δb = b1 - b2

    S_L = asarray(divide(0.040975 * L1, 1 + 0.01765 * L1, where=L1 >= 16))
    S_L[L1 < 16] = 0.511

    S_C = divide(0.0638 * C1, 1 + 0.0131 * C1)

    H1 = arctan2(b1, a1) % (2 * pi)

    T = asarray(0.56 + abs(0.2 * cos(H1 + 0.93 * pi)))
    T[H1 < 0.9 * pi] = 0.36 + abs(0.4 * cos(H1 + 0.19 * pi))
    T[H1 > 1.9 * pi] = 0.36 + abs(0.4 * cos(H1 + 0.19 * pi))

    F = sqrt(divide(C1**4, C1**4 + 1900))

    S_H = S_C * (F * T + 1 - F)

    ΔH = sqrt(Δa**2 + Δb**2 - ΔC**2)

    return sqrt((ΔL / (l * S_L)) ** 2 + (ΔC / (c * S_C)) ** 2 + (ΔH / S_H) ** 2)


def delta_e_cie2000_from_Lab(Lab1, Lab2):
    L1, a1, b1 = moveaxis(Lab1, -1, 0)
    L2, a2, b2 = moveaxis(Lab2, -1, 0)

    C1 = LCh_from_Lab(Lab1)[..., 1]
    C2 = LCh_from_Lab(Lab2)[..., 1]

    C_ = (C1 + C2) / 2

    G = (3 - sqrt(C_**7 / (C_**7 + 25**7))) / 2
    a1 *= G
    a2 *= G

    LCh1 = LCh_from_Lab(moveaxis([L1, a1, b1], 0, -1))
    LCh2 = LCh_from_Lab(moveaxis([L2, a2, b2], 0, -1))

    return delta_e_cie2000_from_LCh(LCh1, LCh2)


def delta_e_cie2000_from_LCh(LCh1, LCh2):
    """Calculates CIEDE2000 color distance between two CIE L*a*b* colors"""
    "https://github.com/lovro-i/CIEDE2000/blob/master/ciede2000.py"

    L1, C1, H1 = moveaxis(LCh1, -1, 0)
    L2, C2, H2 = moveaxis(LCh2, -1, 0)

    L_ = (L1 + L2) / 2
    ΔL = L2 - L1

    C_ = (C1 + C2) / 2
    ΔC = C2 - C1

    H_ = asarray((H1 + H2) / 2)
    Δh = asarray(H2 - H1)

    T = (
        1
        - 0.17 * cos(H_ - pi / 6)
        + 0.24 * cos(2 * H_)
        + 0.32 * cos(3 * H_ + pi / 30)
        - 0.20 * cos(4 * H_ - pi / 1.1)
    )
    # 1.1 ≈ 1.0995574287564276

    Δh[abs(Δh) > pi] %= 2 * pi
    ΔH = 2 * sqrt(C1 * C2) * sin(Δh / 2)

    S_L = 1 + 0.015 * (L_ - 50) ** 2 / sqrt(20 + (L_ - 50) ** 2)
    S_C = 1 + 0.045 * C_
    S_H = 1 + 0.015 * T

    Δθ = exp(-(((H_ - 2.4) / 25) ** 2)) * pi / 6

    R_C = 2 * sqrt(C_**7 / (C_**7 + 25**7))
    R_T = -R_C * sin(2 * Δθ)

    K_L = 1
    K_C = 1
    K_H = 1

    return sqrt(
        (ΔL / (K_L * S_L)) ** 2
        + (ΔC / (K_C * S_C)) ** 2
        + (ΔH / (K_H * S_H)) ** 2
        + R_T * (ΔC / (K_C * S_C)) * (ΔH / (K_H * S_H))
    )
