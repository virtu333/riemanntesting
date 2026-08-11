"""Mod-13 five-wave simplex: even primitive characters, shared gamma factor.

(Z/13)* is cyclic, generator 2.  chi_j(2) = e^{2 pi i j / 12}; parity
chi_j(-1) = (-1)^j, so the even primitive characters are j = 2, 4, 6, 8, 10:
the real quadratic chi_6 and the conjugate pairs (chi_2, chi_10) [order 6]
and (chi_4, chi_8) [order 3].  All five share the gamma factor
(13/pi)^{s/2} Gamma(s/2), so their rotated completed functions

    R_j(s) = e^{-i delta_j} (13/pi)^{s/2} Gamma(s/2) L(s, chi_j),
    delta_j = arg(tau(chi_j)/sqrt(13)) / 2,

all satisfy the SAME Schwarz FE R(1-s) = conj(R(conj(s))) (family lemma,
docs/family_lemma.md, verbatim with the even gamma factor), and every real
combination G = sum c_j W_j of the five real waves W_j(t) = R_j(1/2 + it)
is real on the critical line with mirror-symmetric zeros.

Wave order used everywhere: [chi_2, chi_10, chi_4, chi_8, chi_6]
(indices 0..4; (0,1) and (2,3) are conjugate pairs, 4 is self-dual real).

All five L's are combinations of the SAME 12 Hurwitz zetas zeta(s, a/13):
one batch serves all waves.  Truncation N = 0.55 t + 30 as validated.

Usage: python3 mod13.py     (validation: reality, FE, vertex windings)
"""
import numpy as np
import scipy.special as sp
from dh_core import hurwitz, gamma_np

LOG13PI = np.log(13 / np.pi)
NFAC = 0.55

# index table: 2^k mod 13
_IDX = {}
v = 1
for k in range(12):
    _IDX[v] = k
    v = (v * 2) % 13

EVEN_J = [2, 10, 4, 8, 6]
CHI = {j: {a: np.exp(2j * np.pi * j * _IDX[a] / 12) for a in range(1, 13)}
       for j in EVEN_J}
TAU = {j: sum(CHI[j][a] * np.exp(2j * np.pi * a / 13) for a in range(1, 13))
       for j in EVEN_J}
EPS = {j: TAU[j] / np.sqrt(13) for j in EVEN_J}
DELTA13 = {j: 0.5 * np.angle(EPS[j]) for j in EVEN_J}


def L_all(s, nfac=NFAC):
    """All five even L(s, chi_j) from one Hurwitz batch. Returns list."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    N = int(max(30, nfac * np.max(np.abs(s.imag)) + 30))
    z = {a: hurwitz(s, a / 13.0, N=N) for a in range(1, 13)}
    p = 13.0 ** (-s)
    return [p * sum(CHI[j][a] * z[a] for a in range(1, 13)) for j in EVEN_J]


def theta13(t):
    s = 0.5 + 1j * np.asarray(t, dtype=float)
    return sp.loggamma(s / 2).imag + (np.asarray(t) / 2) * LOG13PI


def waves(ts, chunk=600):
    """The five real waves W_j(t) on the critical line. Returns (5, n)."""
    ts = np.atleast_1d(np.asarray(ts, dtype=float))
    out = np.empty((5, len(ts)))
    for i in range(0, len(ts), chunk):
        t = ts[i:i + chunk]
        Ls = L_all(0.5 + 1j * t)
        th = theta13(t)
        for k, j in enumerate(EVEN_J):
            out[k, i:i + chunk] = (np.exp(1j * (th - DELTA13[j])) * Ls[k]).real
    return out


def pair_factory(i, j):
    """A Z_pair-compatible evaluator returning waves i and j."""
    def pair(ts, chunk=600):
        w = waves(ts, chunk)
        return w[i], w[j]
    return pair


# ---- winding machinery for real combinations (even gamma factor) ----

def f_combo(s, coef):
    Ls = L_all(s)
    return sum(c * np.exp(-1j * DELTA13[j]) * L
               for c, j, L in zip(coef, EVEN_J, Ls) if c != 0)


def phase_F13(s, coef):
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    out = np.empty(len(s))
    for i in range(0, len(s), 600):
        ss = s[i:i + 600]
        lg = sp.loggamma(ss / 2)
        out[i:i + 600] = (lg.imag + (ss / 2).imag * LOG13PI
                          + np.angle(f_combo(ss, coef)))
    return out


def rect_pts(s1, s2, t1, t2, dv, dh_):
    e1 = complex(s1, t1) + np.arange(0, s2 - s1, dh_)
    e2 = complex(s2, t1) + 1j * np.arange(0, t2 - t1, dv)
    e3 = complex(s2, t2) - np.arange(0, s2 - s1, dh_)
    e4 = complex(s1, t2) - 1j * np.arange(0, t2 - t1, dv)
    return np.concatenate([e1, e2, e3, e4, [complex(s1, t1)]])


def _wind(pts, phases):
    d = np.diff(phases)
    d = (d + np.pi) % (2 * np.pi) - np.pi
    return d.sum() / (2 * np.pi), np.max(np.abs(d))


def count_strip13(coef, t1, t2, s1=-1.5, s2=2.5, dv=0.05, dh_=0.02):
    for _ in range(4):
        pts = rect_pts(s1, s2, t1, t2, dv, dh_)
        w, mx = _wind(pts, phase_F13(pts, coef))
        if mx < 1.5 and abs(w - round(w)) < 0.15:
            return int(round(w))
        dv /= 3; dh_ /= 3
    return int(round(w))


def count_box13(coef, s1, s2, t1, t2, dv=0.04, dh_=0.02):
    for _ in range(4):
        pts = rect_pts(s1, s2, t1, t2, dv, dh_)
        w, mx = _wind(pts, np.angle(f_combo(pts, coef)))
        if mx < 1.5 and abs(w - round(w)) < 0.15:
            return int(round(w))
        dv /= 3; dh_ /= 3
    return int(round(w))


def localize_right13(coef, t1, t2):
    found = []

    def rec(a, b, sl, sr, n):
        if n == 0:
            return
        if (b - a) < 0.02 and (sr - sl) < 0.02:
            z = complex(0.5 * (sl + sr), 0.5 * (a + b))
            for _ in range(40):
                h = 1e-6
                fz = f_combo(np.array([z]), coef)[0]
                dfz = (f_combo(np.array([z + h]), coef)[0]
                       - f_combo(np.array([z - h]), coef)[0]) / (2 * h)
                step = fz / dfz
                z -= step
                if abs(step) < 1e-10:
                    break
            found.append((z.real, z.imag))
            return
        if (b - a) >= (sr - sl):
            m = 0.5 * (a + b)
            n1 = count_box13(coef, sl, sr, a, m)
            rec(a, m, sl, sr, n1); rec(m, b, sl, sr, n - n1)
        else:
            m = 0.5 * (sl + sr)
            n1 = count_box13(coef, sl, m, a, b)
            rec(a, b, sl, m, n1); rec(a, b, m, sr, n - n1)

    n = count_box13(coef, 0.501, 2.5, t1, t2)
    rec(t1, t2, 0.501, 2.5, n)
    return n, found


if __name__ == "__main__":
    rng = np.random.default_rng(3)
    print("deltas:", {j: round(DELTA13[j], 6) for j in EVEN_J})
    print("|eps|:", {j: round(abs(EPS[j]), 12) for j in EVEN_J})
    # reality of the five waves
    ts = rng.uniform(10, 300, 5)
    Ls = L_all(0.5 + 1j * ts)
    th = theta13(ts)
    for k, j in enumerate(EVEN_J):
        w = np.exp(1j * (th - DELTA13[j])) * Ls[k]
        print(f"wave chi_{j}: max |Im/Re| = "
              f"{np.max(np.abs(w.imag / w.real)):.2e}")
    # Schwarz FE for random real combos at random strip points
    def R(s, coef):
        return sum(c * np.exp(-1j * DELTA13[j])
                   * (13 / np.pi) ** (s / 2) * gamma_np(s / 2) * L
                   for c, j, L in zip(coef, EVEN_J, L_all(s)))
    ss = rng.uniform(-0.4, 1.4, 6) + 1j * rng.uniform(5, 200, 6)
    for coef in ([1, 1, 0, 0, 0], [1, 0.7, 0.3, 0.2, 0.5],
                 [0, 0, 1, 1, 0], [0.2, 0.2, 0.4, 0.4, 1.0]):
        lhs = R(1 - ss, coef)
        rhs = np.conj(R(np.conj(ss), coef))
        print(f"combo {coef}: FE resid "
              f"{np.max(np.abs(lhs - rhs) / np.abs(lhs)):.2e}")
    # vertices: no right-half zeros for each single L in a sample window
    for k, j in enumerate(EVEN_J):
        coef = [0] * 5
        coef[k] = 1
        n = count_box13(coef, 0.501, 2.5, 100.0, 112.5)
        print(f"vertex chi_{j}: right-half zeros in [100,112.5]: {n}")
    # bookkeeping for a generic combo: strip total vs line count
    coef = [1, 1, 1, 1, 1]
    t1, t2 = 100.0, 112.5
    ts = np.arange(t1, t2, 0.01)
    Z = waves(ts).sum(axis=0)
    nline = int(np.sum(np.signbit(Z[:-1]) != np.signbit(Z[1:])))
    ntot = count_strip13(coef, t1, t2)
    print(f"combo [1,1,1,1,1] in [{t1},{t2}]: line={nline} total={ntot} "
          f"off={ntot - nline}")
