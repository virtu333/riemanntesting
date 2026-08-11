"""Numerical verification of docs/family_lemma.md.

1. Lemma 1: the rotated constituents are real on the critical line.
2. Lemma 2: F_r(1-s) = conj(F_r(conj(s))) for real r, random s in the strip.
3. r = 0 vertex: right-half winding of L_chi alone is 0 in sample windows
   (no detected off-line zeros for the Euler-product vertex — GRH itself is
   of course not checked, only the verified range).

Usage: python3 family.py    (all residuals must print < 1e-10)
"""
import numpy as np
import scipy.special as sp
from dh_core import L_chi, DELTA, gamma_np

LOG5PI = np.log(5 / np.pi)


def Lambda(s, conj=False):
    return (5 / np.pi) ** ((s + 1) / 2) * gamma_np((s + 1) / 2) \
        * L_chi(s, conj=conj)


def F_r(s, r):
    return (np.exp(-1j * DELTA) * Lambda(s)
            + r * np.exp(1j * DELTA) * Lambda(s, conj=True))


def main():
    rng = np.random.default_rng(7)
    # Lemma 1: reality on the line
    ts = rng.uniform(10, 400, 6)
    z1 = np.exp(-1j * DELTA) * Lambda(0.5 + 1j * ts)
    z2 = np.exp(1j * DELTA) * Lambda(0.5 + 1j * ts, conj=True)
    print("Lemma 1  max |Im/Re| Ztilde_1:", np.max(np.abs(z1.imag / z1.real)))
    print("Lemma 1  max |Im/Re| Ztilde_2:", np.max(np.abs(z2.imag / z2.real)))
    # Lemma 2: Schwarz FE for several real r at random strip points
    ss = rng.uniform(-0.4, 1.4, 8) + 1j * rng.uniform(5, 300, 8)
    for r in (0.3, 1.0, 2.5):
        lhs = F_r(1 - ss, r)
        rhs = np.conj(F_r(np.conj(ss), r))
        resid = np.max(np.abs(lhs - rhs) / np.abs(lhs))
        print(f"Lemma 2  r={r}: max FE rel residual = {resid:.2e}")
    # r = 0 vertex: right-half winding of L_chi alone over sample windows
    import hunt2

    def wind_Lchi(s1, s2, t1, t2, dv=0.04, dh_=0.02):
        for _ in range(4):
            pts = hunt2.rect_pts(s1, s2, t1, t2, dv, dh_)
            w, mx = hunt2.wind_closed(pts, np.angle(L_chi(pts)))
            if mx < 1.5 and abs(w - round(w)) < 0.15:
                return int(round(w))
            dv /= 3; dh_ /= 3
        return int(round(w))

    for (t1, t2) in ((80.0, 92.5), (300.0, 312.5), (480.0, 492.5)):
        n = wind_Lchi(0.501, 2.5, t1, t2)
        print(f"vertex r=0: right-half zeros of L_chi in [{t1},{t2}]: {n}")


if __name__ == "__main__":
    main()
