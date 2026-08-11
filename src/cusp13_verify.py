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
    bZ2, bZ3 = float(z[1, 0]), float(z[2, 0])       # constant-term coeffs
    # tangent (fold) direction kills the constant term: drho = -(Z2/Z3) dr;
    # along it the linear term has coefficient gamma_eff = -W(Z2, Z3)/Z3
    gamma_eff = float((d1[1, 0] * z[2, 0] - z[1, 0] * d1[2, 0]) / z[2, 0])
    print(f"alpha=G'''={alpha:.4f}  Z2(tc)={bZ2:.4f}  Z3(tc)={bZ3:.4f}  "
          f"gamma_eff={gamma_eff:.4f}")
    w = 0.5 * np.log(13 * tc / (2 * np.pi))
    msp = np.pi / w
    t1, t2 = tc - 3 * msp, tc + 3 * msp
    # at the cusp: triple zero = 1 sign change but 3 winding counts -> off=2
    off0, _ = probe([1.0, rc, rho], t1, t2)
    print(f"at the cusp: off={off0} (expect 2: triple zero counts 3 in the "
          f"strip, 1 as a sign change — multiplicity, not off-line zeros)")
    rays = (
        ("rho (generic)", lambda dd: [1.0, rc, rho + dd],
         lambda dd: abs(bZ3 * dd), 1 / 3),
        ("r (generic)", lambda dd: [1.0, rc + dd, rho],
         lambda dd: abs(bZ2 * dd), 1 / 3),
        ("tangent (fold ctrl)",
         lambda dd: [1.0, rc + dd, rho - dd * bZ2 / bZ3],
         None, 1 / 2),
    )
    for name, mk, c0, expo in rays:
        print(f"\n{name}-ray, expected exponent {expo:.3f}:")
        rows = []
        for dd in (0.003, 0.01, 0.03, 0.1):
            for sgn in (+1, -1):
                d = sgn * dd
                if c0 is not None:
                    pred = (np.sqrt(3) / 2) * abs(6 * c0(d) / alpha) ** (1 / 3)
                else:
                    pred = abs(6 * gamma_eff * d / alpha) ** 0.5 \
                        if gamma_eff * d / alpha > 0 else 0.0
                off, zs = probe(mk(d), t1, t2)
                exc = max((s - 0.5 for s, _ in zs), default=0.0)
                rows.append((dd, off, exc))
                print(f"  d={d:+7.3f}: off={off} exc={exc:.4f} "
                      f"pred~{pred:.4f}", flush=True)
        pts = [(np.log(dd), np.log(exc)) for dd, off, exc in rows if exc > 0]
        if len(pts) >= 3:
            x, y = np.array(pts).T
            print(f"  fitted exponent: {np.polyfit(x, y, 1)[0]:.3f}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
