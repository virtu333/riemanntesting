"""Fast dual line-scan of both DH constituents (Step 2 workhorse).

Both L(s,chi) and L(s,chibar) are linear combinations of the SAME four Hurwitz
zetas zeta(s, a/5), so one batch evaluation serves both waves (2x saving), and
Euler-Maclaurin truncation N = 0.55 t keeps ~1e-12 relative accuracy on the
critical line (validated against mpmath at t = 1000/3000/5000; the 1.1 t
default in dh_core is very conservative) for another ~2.5x.

Per window: grid scan at DT, sign changes -> vectorized bisection (both waves'
brackets in one batch), slopes by central difference.  Missed-zero insurance:
 1. slope-alternation check per wave (consecutive zeros of a continuous
    function must alternate slope signs; a break = one missed zero) -> local
    rescan at DT/20;
 2. dip detector: grid points where |Z| dips below 0.02 x local scale with no
    sign change within +-3 grid steps could hide a near-tangent PAIR of zeros
    (alternation-safe miss) -> same local rescan.

Usage: python3 scan2.py T1 T2 TAG    (writes ../data/scan_<TAG>.npz with
                                      g1, d1, g2, d2 for zeros in [T1, T2))
"""
import numpy as np
import scipy.special as sp
import sys
from dh_core import CHI, DELTA, hurwitz

LOG5PI = np.log(5 / np.pi)
NFAC = 0.55
DT = 0.02
H = 1e-4


def L_pair(s, nfac=NFAC):
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    N = int(max(30, nfac * np.max(np.abs(s.imag)) + 30))
    z = {a: hurwitz(s, a / 5.0, N=N) for a in (1, 2, 3, 4)}
    L1 = sum(CHI[a] * z[a] for a in (1, 2, 3, 4))
    L2 = sum(np.conj(CHI[a]) * z[a] for a in (1, 2, 3, 4))
    p = 5.0 ** (-s)
    return p * L1, p * L2


def Z_pair(ts, chunk=800):
    """Both real rotated Z-functions on the critical line, one L evaluation."""
    ts = np.atleast_1d(np.asarray(ts, dtype=float))
    Z1, Z2 = np.empty(len(ts)), np.empty(len(ts))
    for i in range(0, len(ts), chunk):
        t = ts[i:i + chunk]
        s = 0.5 + 1j * t
        lg = sp.loggamma((s + 1) / 2)
        base = lg.imag + (t / 2) * LOG5PI
        L1, L2 = L_pair(s)
        Z1[i:i + chunk] = (np.exp(1j * (base - DELTA)) * L1).real
        Z2[i:i + chunk] = (np.exp(1j * (base + DELTA)) * L2).real
    return Z1, Z2


def bisect_batch(a, b, fa, wave, iters=30):
    """Vectorized bisection on one wave's brackets."""
    a, b, fa = a.copy(), b.copy(), fa.copy()
    for _ in range(iters):
        m = 0.5 * (a + b)
        Zm = Z_pair(m)[wave]
        left = np.signbit(Zm) == np.signbit(fa)
        a = np.where(left, m, a)
        fa = np.where(left, Zm, fa)
        b = np.where(left, b, m)
    return 0.5 * (a + b)


def scan_window(t1, t2, dt=DT):
    """Zeros and slopes of both waves in [t1, t2). Grid overshoots t2 by one
    spacing so brackets straddling t2 are caught by exactly one window."""
    ts = np.arange(t1, t2 + dt, dt)
    Z1, Z2 = Z_pair(ts)
    out = []
    for wave, Z in ((0, Z1), (1, Z2)):
        idx = np.where(np.signbit(Z[:-1]) != np.signbit(Z[1:]))[0]
        idx = idx[ts[idx] < t2]
        zs = bisect_batch(ts[idx], ts[idx + 1], Z[idx], wave)
        # dip detector: deep |Z| minima with no sign change nearby
        near = np.zeros(len(ts), dtype=bool)
        for k in idx:
            near[max(0, k - 3):k + 5] = True
        absZ = np.abs(Z)
        interior = np.zeros(len(ts), dtype=bool)
        interior[1:-1] = (absZ[1:-1] < absZ[:-2]) & (absZ[1:-1] < absZ[2:])
        scale = np.maximum(np.max(absZ), 1e-6)
        sus = np.where(interior & ~near & (absZ < 0.02 * scale))[0]
        for k in sus:
            tf = np.arange(ts[k] - 2 * dt, ts[k] + 2 * dt, dt / 20)
            Zf = Z_pair(tf)[wave]
            j = np.where(np.signbit(Zf[:-1]) != np.signbit(Zf[1:]))[0]
            if len(j):
                zs = np.sort(np.concatenate(
                    [zs, bisect_batch(tf[j], tf[j + 1], Zf[j], wave, 20)]))
        if len(zs) > 1:  # rescans can re-find an already-known zero: dedup
            zs = zs[np.concatenate([[True], np.diff(zs) > 1e-7])]
        d = (Z_pair(zs + H)[wave] - Z_pair(zs - H)[wave]) / (2 * H)
        # slope-alternation check: a break means a missed zero -> rescan
        bad = np.where(np.sign(d[1:]) == np.sign(d[:-1]))[0]
        for k in bad:
            tf = np.arange(zs[k], zs[k + 1], dt / 20)
            Zf = Z_pair(tf)[wave]
            j = np.where(np.signbit(Zf[:-1]) != np.signbit(Zf[1:]))[0]
            if len(j):
                extra = bisect_batch(tf[j], tf[j + 1], Zf[j], wave, 20)
                zs = np.sort(np.concatenate([zs, extra]))
                zs = zs[np.concatenate([[True], np.diff(zs) > 1e-7])]
                d = (Z_pair(zs + H)[wave] - Z_pair(zs - H)[wave]) / (2 * H)
        if np.any(np.sign(d[1:]) == np.sign(d[:-1])):
            print(f"  WARNING wave {wave}: alternation break persists "
                  f"in [{t1},{t2})", flush=True)
        out.append((zs, d))
    return out


if __name__ == "__main__":
    T1, T2, tag = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    g1s, d1s, g2s, d2s = [], [], [], []
    for t1 in np.arange(T1, T2, 12.5):
        t2 = min(t1 + 12.5, T2)
        (z1, s1), (z2, s2) = scan_window(t1, t2)
        g1s.append(z1); d1s.append(s1); g2s.append(z2); d2s.append(s2)
        print(f"[{t1:7.1f},{t2:7.1f}] wave1 {len(z1):3d}  wave2 {len(z2):3d}",
              flush=True)
    np.savez(f"../data/scan_{tag}.npz",
             g1=np.concatenate(g1s), d1=np.concatenate(d1s),
             g2=np.concatenate(g2s), d2=np.concatenate(d2s))
    print("saved", f"../data/scan_{tag}.npz")
