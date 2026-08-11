"""Wronskian bifurcation census for the family G_r(t) = Z1(t) + r Z2(t).

A pair of on-line zeros of G_r leaves (or lands on) the critical line
through a double zero: G_r(t*) = G_r'(t*) = 0.  Eliminating r gives the
scalar condition on the two fixed constituents

    H(t) = Z1(t) Z2'(t) - Z1'(t) Z2(t) = 0,

and every simple zero t* of H yields the critical coefficient

    r* = -Z1(t*)/Z2(t*)  ( = -Z1'(t*)/Z2'(t*), used as a consistency check).

In the two-wave picture H ~ A B omega sin(dphi): its zeros at dphi = pi are
the anti-phase alignments (r* > 0, our family), at dphi = 0 the in-phase
alignments (r* < 0, the mirror family Z1 - |r| Z2) — the Step-1 necessity
law is the sign structure of the Wronskian.

Fold normal form at a census point (Taylor of G_r in both variables):

    G_r(t) ~ Z2(t*) (r - r*) + 1/2 G''_{r*}(t*) (t - t*)^2
    =>  (t - t*)^2 ~ -2 Z2(t*) / G''(t*) * (r - r*),

so the pair is on-line for one sign of (r - r*) and a complex-conjugate
pair (off-line, |sigma - 1/2| = |Im t|) for the other:

    |sigma - 1/2| = |2 Z2(t*) / G''(t*)|^{1/2} |r - r*|^{1/2}.

Per census point we record: t*, r* (both formulas), the fold coefficient
C = |2 Z2 / G''|^{1/2}, the side of r* on which the pair is complex
(sign of Z2/G''; complex side is where Z2/G'' * (r - r*) > 0 ... i.e.
(t-t*)^2 < 0), and the degeneracy diagnostic |G''| / ((|A|+|B|) omega^2)
(near-zero flags a cluster-mediated / almost-cusp event: a genuine cusp
G = G' = G'' = 0 is nongeneric in a one-parameter family and is NOT claimed
— see review notes in FINDINGS).

Method: dense grid (DT = 0.02) of Z1, Z2 via the dual evaluator; H on the
grid with grid-spline derivatives to locate sign changes; each root refined
by bisection on H evaluated with accurate central differences (h = 1e-4);
derivatives of G at the root by 5-point stencils.

Usage: python3 bifurcate.py T1 T2 TAG    (writes ../data/census_<TAG>.npz)
"""
import sys
import numpy as np
from scan2 import Z_pair

DT = 0.02
H_STEP = 1e-4     # stencil for accurate Z' at refinement stage


def H_accurate(ts):
    """H = Z1 Z2' - Z1' Z2 with central-difference derivatives, batched."""
    ts = np.atleast_1d(np.asarray(ts, dtype=float))
    Zp = Z_pair(np.concatenate([ts, ts + H_STEP, ts - H_STEP]))
    n = len(ts)
    z1, z2 = Zp[0][:n], Zp[1][:n]
    d1 = (Zp[0][n:2 * n] - Zp[0][2 * n:]) / (2 * H_STEP)
    d2 = (Zp[1][n:2 * n] - Zp[1][2 * n:]) / (2 * H_STEP)
    return z1 * d2 - d1 * z2, z1, z2, d1, d2


def census_window(t1, t2):
    ts = np.arange(t1, t2 + DT, DT)
    Z1, Z2 = Z_pair(ts)
    # grid derivatives (spline-quality is enough to locate sign changes;
    # roots are refined with accurate stencils afterwards)
    d1 = np.gradient(Z1, DT)
    d2 = np.gradient(Z2, DT)
    Hg = Z1 * d2 - d1 * Z2
    idx = np.where(np.signbit(Hg[:-1]) != np.signbit(Hg[1:]))[0]
    idx = idx[ts[idx] < t2]
    if not len(idx):
        return []
    # vectorized bisection on accurate H
    a, b = ts[idx].copy(), ts[idx + 1].copy()
    fa = H_accurate(a)[0]
    for _ in range(40):
        m = 0.5 * (a + b)
        fm = H_accurate(m)[0]
        left = np.signbit(fm) == np.signbit(fa)
        a = np.where(left, m, a)
        fa = np.where(left, fm, fa)
        b = np.where(left, b, m)
    tstar = 0.5 * (a + b)
    # values and derivatives at the roots: 5-point stencil for G''
    h = 1e-3
    grid = np.concatenate([tstar + k * h for k in (-2, -1, 0, 1, 2)])
    Z1g, Z2g = Z_pair(grid)
    n = len(tstar)
    z1s = {k: Z1g[(k + 2) * n:(k + 3) * n] for k in (-2, -1, 0, 1, 2)}
    z2s = {k: Z2g[(k + 2) * n:(k + 3) * n] for k in (-2, -1, 0, 1, 2)}
    out = []
    for i in range(n):
        z1, z2 = z1s[0][i], z2s[0][i]
        d1_ = (z1s[1][i] - z1s[-1][i]) / (2 * h)
        d2_ = (z2s[1][i] - z2s[-1][i]) / (2 * h)
        dd1 = (z1s[1][i] - 2 * z1 + z1s[-1][i]) / h ** 2
        dd2 = (z2s[1][i] - 2 * z2 + z2s[-1][i]) / h ** 2
        if abs(z2) < 1e-12 or abs(d2_) < 1e-12:
            continue  # degenerate cross-collision: both formulas blow up
        r1, r2 = -z1 / z2, -d1_ / d2_
        rstar = r1
        gpp = dd1 + rstar * dd2
        w = 0.5 * np.log(5 * tstar[i] / (2 * np.pi))
        # local amplitude scale from the oscillator envelopes
        A = np.hypot(z1, d1_ / w)
        B = np.hypot(z2, d2_ / w)
        out.append(dict(
            t=float(tstar[i]), r=float(rstar), r_check=float(r2),
            rerr=float(abs(r1 - r2) / max(abs(r1), abs(r2), 1e-30)),
            C=float(np.sqrt(abs(2 * z2 / gpp))) if gpp != 0 else np.inf,
            # complex side: (t-t*)^2 = -2 z2/gpp (r - r*) < 0
            #   => complex for sign(r - r*) = sign(z2/gpp)
            side=int(np.sign(z2 / gpp)) if gpp != 0 else 0,
            degen=float(abs(gpp) / ((A + B) * w ** 2)),
            A=float(A), B=float(B), omega=float(w),
        ))
    return out


if __name__ == "__main__":
    T1, T2, tag = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    rows = []
    for t1 in np.arange(T1, T2, 12.5):
        t2 = min(t1 + 12.5, T2)
        rows += census_window(t1, t2)
        print(f"[{t1:7.1f},{t2:7.1f}] census points so far: {len(rows)}",
              flush=True)
    np.savez(f"../data/census_{tag}.npz",
             **{k: np.array([r[k] for r in rows])
                for k in rows[0].keys()})
    npos = sum(1 for r in rows if 0 < r["r"] <= 1)
    print(f"saved ../data/census_{tag}.npz: {len(rows)} points, "
          f"{npos} with 0 < r* <= 1")
