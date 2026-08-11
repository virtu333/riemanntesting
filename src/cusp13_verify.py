"""Verify a census cusp: cube-root vs square-root unfolding scaling.

At a cusp (t_c, r_c, rho_c) of G = Z1 + r Z2 + rho Z3 the normal form is

    G ~ (alpha/6) tau^3 + gamma dr tau + beta drho,
    tau = t - t_c,  alpha = G'''(t_c),  beta = Z3(t_c),  gamma = Z2'(t_c),

(the tau^2 term vanishes by G'' = 0; Z2(t_c) contributes to beta-like and
Z3'(t_c) to gamma-like terms at higher order in the small parameters).
Predictions tested by strip winding of the exact family:
 - rho-ray (dr = 0): one real zero + a complex pair BOTH sides, excursion
       |Im t| = (sqrt(3)/2) |6 beta drho / alpha|^{1/3}   (CUBE ROOT);
 - r-ray (drho = 0): fold-pair behavior, three real zeros on one side,
   one real + complex pair on the other with
       |Im t| ~ |2 gamma dr / alpha|^{1/2} * sqrt(3)/..   (SQUARE ROOT;
   slope is the tested content, the prefactor for this ray is quoted to
   O(1)).
The slope of log(excursion) vs log(perturbation) is the universality-class
signature: 1/3 on the rho-ray vs 1/2 on the r-ray.

Usage: python3 cusp13_verify.py [k]   (k-th best cusp; default: the cusp
       with smallest svd residual and |r_c|, |rho_c| in [0.05, 20])
"""
import sys
import numpy as np
import mod13
from cusp13 import stencil, IDX


def line_count13(coef, t1, t2, dt=0.01):
    ts = np.arange(t1, t2, dt)
    W = mod13.waves(ts)
    Z = coef[0] * W[IDX[0]] + coef[1] * W[IDX[1]] + coef[2] * W[IDX[2]]
    return int(np.sum(np.signbit(Z[:-1]) != np.signbit(Z[1:]))), ts, Z


def full_coef(c3):
    coef = [0.0] * 5
    for k, idx in enumerate(IDX):
        coef[idx] = c3[k]
    return coef


def probe(c3, t1, t2):
    coef = full_coef(c3)
    nline, ts, Z = line_count13(c3, t1, t2)
    ntot = mod13.count_strip13(coef, t1, t2)
    off = ntot - nline
    zs = []
    if off > 0 and off % 2 == 0:
        nr, zs = mod13.localize_right13(coef, t1, t2)
    return off, zs


def main(pick=None):
    c = np.load("../data/cusps13.npz")
    ok = (c["boundary"] == 0) & (np.abs(c["r"]) > 0.05) \
        & (np.abs(c["r"]) < 20) & (np.abs(c["rho"]) > 0.05) \
        & (np.abs(c["rho"]) < 20) & (c["t"] > 100)
    order = np.argsort(c["sing"][ok])
    idx = np.where(ok)[0][order[pick or 0]]
    tc, rc, rho = c["t"][idx], c["r"][idx], c["rho"][idx]
    print(f"cusp: t_c={tc:.6f} r_c={rc:.6f} rho_c={rho:.6f} "
          f"sing={c['sing'][idx]:.2e} g3={c['g3'][idx]:.4f}")
    z, d1, _, d3 = stencil(np.array([tc]))
    alpha = float(d3[0, 0] + rc * d3[1, 0] + rho * d3[2, 0])
    beta = float(z[2, 0])
    gamma = float(d1[1, 0])
    print(f"alpha=G'''={alpha:.4f}  beta=Z3={beta:.4f}  gamma=Z2'={gamma:.4f}")
    w = 0.5 * np.log(13 * tc / (2 * np.pi))
    msp = np.pi / w
    t1, t2 = tc - 3 * msp, tc + 3 * msp
    # sanity at the cusp point itself
    off0, _ = probe([1.0, rc, rho], t1, t2)
    print(f"at the cusp: off={off0} (expect 0: the triple zero is ON line)")
    for name, ray, expo in (("rho", "rho", 1 / 3), ("r", "r", 1 / 2)):
        print(f"\n{name}-ray (expected exponent {expo:.3f}):")
        rows = []
        for dd in (0.003, 0.01, 0.03, 0.1):
            for sgn in (+1, -1):
                if ray == "rho":
                    c3 = [1.0, rc, rho + sgn * dd]
                    pred = (np.sqrt(3) / 2) * abs(6 * beta * dd / alpha) ** (1 / 3)
                else:
                    c3 = [1.0, rc + sgn * dd, rho]
                    pred = abs(2 * gamma * dd / alpha) ** 0.5
                off, zs = probe(c3, t1, t2)
                exc = max((s - 0.5 for s, _ in zs), default=0.0)
                rows.append((dd, sgn, off, exc, pred))
                print(f"  d{ray}={sgn*dd:+7.3f}: off={off} "
                      f"exc={exc:.4f} pred~{pred:.4f}", flush=True)
        # slope fit on the side(s) with a complex pair
        pts = [(np.log(dd), np.log(exc)) for dd, sgn, off, exc, _ in rows
               if exc > 0]
        if len(pts) >= 3:
            x, y = np.array(pts).T
            print(f"  fitted exponent: {np.polyfit(x, y, 1)[0]:.3f}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
