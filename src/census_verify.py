"""Direct validation of the fold normal form at census points.

For a census fold (t*, r*, C): at r = r* + delta on the complex side the
family G_r has an off-line pair with |sigma - 1/2| ~ C sqrt(delta); at
r = r* - delta on the real side it has an on-line pair with gap
|t+ - t-| ~ 2 C sqrt(delta).  Both are measured here with the exact family
(strip winding + localization for the complex side, line scan for the real
side) and compared to the prediction.

Usage: python3 census_verify.py     (folds and deltas hardcoded below)
"""
import numpy as np
import scan2
import hunt2
from dh_core import DELTA


def set_r(r):
    """Point hunt2's winding machinery at f_r = e^{-i d} L1 + r e^{i d} L2."""
    def f_r(s):
        L1, L2 = scan2.L_pair(s)
        return np.exp(-1j * DELTA) * L1 + r * np.exp(1j * DELTA) * L2
    hunt2.f_dh = f_r


def G_r(ts, r):
    Z1, Z2 = scan2.Z_pair(ts)
    return Z1 + r * Z2


def real_side_gap(tstar, r):
    ts = np.arange(tstar - 2.0, tstar + 2.0, 0.002)
    Z = G_r(ts, r)
    sc = ts[np.where(np.signbit(Z[:-1]) != np.signbit(Z[1:]))[0]]
    if len(sc) < 2:
        return None
    i = np.argsort(np.abs(sc - tstar))[:2]
    z1, z2 = sorted(sc[i])
    return z2 - z1


def complex_side_zero(tstar, r):
    """Wind [tstar - 2, tstar + 2] with f_r, localize the right-half zero."""
    set_r(r)
    ts = np.arange(tstar - 3.0, tstar + 3.0, 0.01)
    Z = G_r(ts, r)
    zt = ts[np.where(np.signbit(Z[:-1]) != np.signbit(Z[1:]))[0]] + 0.005

    def snap(edge):
        k = np.searchsorted(zt, edge)
        if k == 0 or k == len(zt):
            return edge
        return 0.5 * (zt[k - 1] + zt[k])

    a, b = snap(tstar - 2.0), snap(tstar + 2.0)
    nline = int(np.sum((zt > a) & (zt < b)))
    ntot = hunt2.count_strip(a, b)
    off = ntot - nline
    if off <= 0:
        return off, None
    nr, zs = hunt2.localize_right(a, b, off // 2)
    z = min(zs, key=lambda z: abs(z[1] - tstar)) if zs else None
    return off, z


FOLDS = [  # (t*, r*, C) from the census
    (85.89902533655524, 0.43974, 0.56833),
    (2452.4559, 0.97440, 0.22235),
    (2923.5580, 0.55351, 0.19301),
]
DELTAS = [0.01, 0.03, 0.1]

if __name__ == "__main__":
    for tstar, rstar, C in FOLDS:
        print(f"fold t*={tstar:.4f} r*={rstar:.4f} C={C:.4f}")
        for d in DELTAS:
            pred = C * np.sqrt(d)
            gap = real_side_gap(tstar, rstar - d)
            gr = gap / (2 * pred) if gap else np.nan
            off, z = complex_side_zero(tstar, rstar + d)
            if z:
                exc = z[0] - 0.5
                print(f"  delta={d:5.3f}: pred={pred:.4f}  "
                      f"real-side gap/2pred={gr:.3f}  "
                      f"complex side off={off} sigma-1/2={exc:.4f} "
                      f"ratio={exc/pred:.3f} (t={z[1]:.4f})")
            else:
                print(f"  delta={d:5.3f}: pred={pred:.4f}  "
                      f"real-side gap/2pred={gr:.3f}  "
                      f"complex side off={off} NO ZERO FOUND")
