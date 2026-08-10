"""Collision-mechanism analysis: Z_DH(t) = Z_chi(t) + Z_chibar(t) exactly, so off-line
zeros of DH are studied via near-collisions between the zero sequences of the two
constituent L-functions on the critical line.

Pipeline:
  1. zeros_of(conj): all zeros of the rotated-real Z-function of L(s,chi) / L(s,chibar)
     on the line up to T (sign-change scan + bisection). Saves g1.npy / g2.npy.
  2. Mutual-nearest cross pairs (dedup: keep (i,j) only if each is the other's nearest).
  3. Per pair: slopes at the two zeros -> sgn = sign(slope product)
       sgn > 0  "in-phase"  (aligned crossings)
       sgn < 0  "anti-phase" (opposed crossings; two-wave cancellation possible)
     amplitudes A, B = |slope|/omega, omega(t) = 0.5*log(5t/2pi),
     mean spacing = pi/omega, cancellation depth kappa = residual amplitude /(A+B).
  4. Outcome label: off-line zero (from hunt catalog) within 1.0 of pair midpoint.
     NOTE: window can double-count — two adjacent collisions can map to one off-line
     zero. For clean counts, assign each off-line zero to its nearest collision only.

Established result at T <= 500 (reproduce before extending):
  304 mutual-nearest pairs; 230 in-phase -> 0 off-line; 74 anti-phase -> all 10
  off-line mirror pairs. corr(sigma-1/2, kappa) ~ 0.75 among off-line.

Usage: python3 collisions.py [T]      (default 500; needs data/offline_zeros json)
"""
import numpy as np, json, sys, os
import scipy.special as sp
from dh_core import L_chi, DELTA

LOG5PI = np.log(5 / np.pi)

def Z_single(ts, conj, chunk=1500):
    """Real rotated Z for a single L-function on the critical line.
    Reality: conj(Lambda) = eps_bar * Lambda on the line => e^{-i delta} Lambda real
    (delta = arg(eps)/2 for chi; -delta for chibar)."""
    ts = np.atleast_1d(np.asarray(ts, dtype=float))
    out = np.empty(len(ts))
    d = DELTA if conj else -DELTA
    for i in range(0, len(ts), chunk):
        t = ts[i:i + chunk]; s = 0.5 + 1j * t
        lg = sp.loggamma((s + 1) / 2)
        ph = lg.imag + (t / 2) * LOG5PI + d
        out[i:i + chunk] = (np.exp(1j * ph) * L_chi(s, conj=conj)).real
    return out

def zeros_of(conj, T=500.0, dt=0.005):
    ts = np.arange(0.5, T, dt)
    Z = Z_single(ts, conj)
    idx = np.where(np.signbit(Z[:-1]) != np.signbit(Z[1:]))[0]
    zs = []
    for i in idx:
        a, b = ts[i], ts[i + 1]
        fa = Z[i]
        for _ in range(40):
            m = 0.5 * (a + b)
            fm = Z_single([m], conj)[0]
            if np.signbit(fm) == np.signbit(fa):
                a = m; fa = fm
            else:
                b = m
        zs.append(0.5 * (a + b))
    return np.array(zs)

def omega(t):
    return 0.5 * np.log(5 * t / (2 * np.pi))

def analyze(g1, g2, off_t, off_s):
    h = 1e-4
    d1 = (Z_single(g1 + h, False) - Z_single(g1 - h, False)) / (2 * h)
    d2 = (Z_single(g2 + h, True) - Z_single(g2 - h, True)) / (2 * h)
    # mutual-nearest cross pairs
    j = np.searchsorted(g2, g1); j = np.clip(j, 1, len(g2) - 1)
    uniq = []
    for i, g in enumerate(g1):
        jj = j[i] if abs(g2[j[i]] - g) < abs(g2[j[i] - 1] - g) else j[i] - 1
        k = np.searchsorted(g1, g2[jj]); k = np.clip(k, 1, len(g1) - 1)
        kk = k if abs(g1[k] - g2[jj]) < abs(g1[k - 1] - g2[jj]) else k - 1
        if kk == i:
            uniq.append((i, jj))
    rows = []
    for i, jj in uniq:
        gap = abs(g1[i] - g2[jj]); mid = 0.5 * (g1[i] + g2[jj])
        w = omega(mid); msp = np.pi / w
        A, B = abs(d1[i]) / w, abs(d2[jj]) / w
        sgn = int(np.sign(d1[i] * d2[jj]))
        dphi = w * gap if sgn > 0 else np.pi - w * gap
        if sgn < 0:
            C = np.sqrt(max(A * A + B * B - 2 * A * B * np.cos(np.pi - dphi), 0))
        else:
            C = np.sqrt(A * A + B * B + 2 * A * B * np.cos(dphi))
        kappa = C / (A + B)
        doff = np.min(np.abs(off_t - mid)) if len(off_t) else np.inf
        out = 1 if doff < 1.0 else 0
        sig = float(off_s[np.argmin(np.abs(off_t - mid))] - 0.5) if out else None
        rows.append(dict(mid=float(mid), gap=float(gap), gapn=float(gap / msp),
                         sgn=sgn, A=float(A), B=float(B),
                         rho=float(min(A, B) / max(A, B)), kappa=float(kappa),
                         out=out, sig=sig))
    return rows

if __name__ == "__main__":
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 500.0
    if os.path.exists("g1.npy"):
        g1, g2 = np.load("g1.npy"), np.load("g2.npy")
    else:
        g1 = zeros_of(False, T); g2 = zeros_of(True, T)
        np.save("g1.npy", g1); np.save("g2.npy", g2)
    print(f"L(chi) zeros: {len(g1)}, L(chibar): {len(g2)}")
    cat = json.load(open("../data/offline_zeros_T500.json"))
    off_s = np.array([z["sigma"] for z in cat]); off_t = np.array([z["t"] for z in cat])
    rows = analyze(g1, g2, off_t, off_s)
    json.dump(rows, open("collisions_full.json", "w"))
    n_anti = sum(1 for r in rows if r["sgn"] < 0)
    for s_ in (-1, 1):
        on = sum(1 for r in rows if r["sgn"] == s_ and not r["out"])
        off = sum(1 for r in rows if r["sgn"] == s_ and r["out"])
        print(f"  sign={s_:+d}: on-line {on:3d}   off-line {off:2d}")
