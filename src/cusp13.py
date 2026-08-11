"""Cusp census for a three-wave mod-13 family G = Z1 + r Z2 + rho Z3.

A genuine cusp G = G' = G'' = 0 with (1, r, rho) real requires the 3x3
Wronskian

    D(t) = det [ [Z1, Z2, Z3], [Z1', Z2', Z3'], [Z1'', Z2'', Z3''] ]

to vanish; (1, r_c, rho_c) is the null vector.  D is a smooth scalar with
no poles (unlike the Cramer parameterization r(t), rho(t), which blows up
on the 2x2-minor zero set), so the census is again a 1-D sign-change scan.
In a one-parameter family cusps are non-generic (verified empirically at
mod 5: no near-degenerate folds); with two coefficients they are isolated
points — this is the instrument for the cusp-exponent question.

Local normal form at a cusp: G ~ (alpha/6)(t-t_c)^3 + beta drho + ..., so
perturbing rho alone gives one real zero plus a complex pair with
|Im t| = (sqrt(3)/2) |6 beta drho / alpha|^{1/3}: CUBE-ROOT scaling,
against the fold's square root.  The scaling test is run separately
(cusp13_verify.py) by winding along a rho-ray.

Family here: waves (0, 1, 4) = [chi_2, chi_10, chi_6]: contains the
DH-analog subfamily at rho = 0 and the real-character vertex.

Usage: python3 cusp13.py T1 T2      (writes ../data/cusps13.npz)
"""
import sys
import numpy as np
import mod13

IDX = (0, 1, 4)
DT = 0.02
H = 1e-3


def stencil(ts):
    """Z, Z', Z'', Z''' for the three waves at ts via 5-point stencils."""
    ts = np.atleast_1d(np.asarray(ts, dtype=float))
    grid = np.concatenate([ts + k * H for k in (-2, -1, 0, 1, 2)])
    W = mod13.waves(grid)[list(IDX), :]
    n = len(ts)
    f = {k: W[:, (k + 2) * n:(k + 3) * n] for k in (-2, -1, 0, 1, 2)}
    z = f[0]
    d1 = (f[-2] - 8 * f[-1] + 8 * f[1] - f[2]) / (12 * H)
    d2 = (f[1] - 2 * f[0] + f[-1]) / H ** 2
    d3 = (-0.5 * f[-2] + f[-1] - f[1] + 0.5 * f[2]) / H ** 3
    return z, d1, d2, d3


def D_accurate(ts):
    z, d1, d2, _ = stencil(ts)
    # z etc have shape (3, n): build (n, 3, 3) with rows = derivative order
    A = np.stack([z.T, d1.T, d2.T], axis=1)    # (n, 3, 3): row=deriv, col=wave
    return np.linalg.det(A)


def main(T1, T2):
    rows = []
    for t1 in np.arange(T1, T2, 12.5):
        t2 = min(t1 + 12.5, T2)
        ts = np.arange(t1, t2 + DT, DT)
        W = mod13.waves(ts)[list(IDX), :]
        d1 = np.gradient(W, DT, axis=1)
        d2 = np.gradient(d1, DT, axis=1)
        A = np.stack([W.T, d1.T, d2.T], axis=1)
        D = np.linalg.det(A)
        idx = np.where(np.signbit(D[:-1]) != np.signbit(D[1:]))[0]
        idx = idx[ts[idx] < t2]
        if len(idx):
            a, b = ts[idx].copy(), ts[idx + 1].copy()
            fa = D_accurate(a)
            for _ in range(35):
                m = 0.5 * (a + b)
                fm = D_accurate(m)
                left = np.signbit(fm) == np.signbit(fa)
                a = np.where(left, m, a)
                fa = np.where(left, fm, fa)
                b = np.where(left, b, m)
            tc = 0.5 * (a + b)
            z, dz1, dz2, dz3 = stencil(tc)
            for i in range(len(tc)):
                M = np.array([[z[0, i], z[1, i], z[2, i]],
                              [dz1[0, i], dz1[1, i], dz1[2, i]],
                              [dz2[0, i], dz2[1, i], dz2[2, i]]])
                _, s, vt = np.linalg.svd(M)
                v = vt[-1]
                if abs(v[0]) < 1e-6 * np.max(np.abs(v)):
                    rows.append(dict(t=float(tc[i]), r=np.nan, rho=np.nan,
                                     sing=float(s[-1] / s[0]), boundary=1))
                    continue
                v = v / v[0]
                rc, rho = float(v[1]), float(v[2])
                g3 = float(dz3[0, i] + rc * dz3[1, i] + rho * dz3[2, i])
                w = 0.5 * np.log(13 * tc[i] / (2 * np.pi))
                amp = sum(np.hypot(z[k, i], dz1[k, i] / w) for k in range(3))
                rows.append(dict(
                    t=float(tc[i]), r=rc, rho=rho,
                    sing=float(s[-1] / s[0]),
                    g3=g3, degen3=float(abs(g3) / (amp * w ** 3)),
                    boundary=0))
        print(f"[{t1:7.1f},{t2:7.1f}] cusp candidates so far: {len(rows)}",
              flush=True)
    np.savez("../data/cusps13.npz",
             **{k: np.array([r.get(k, np.nan) for r in rows])
                for k in ("t", "r", "rho", "sing", "g3", "degen3", "boundary")})
    good = [r for r in rows if not r["boundary"]]
    both_pos = [r for r in good if r["r"] > 0 and r["rho"] > 0]
    print(f"saved ../data/cusps13.npz: {len(rows)} candidates, "
          f"{len(good)} with valid (r, rho), {len(both_pos)} in the "
          f"positive quadrant")
    for r in sorted(good, key=lambda x: x["t"])[:15]:
        print(f"  t_c={r['t']:9.3f} r_c={r['r']:+9.4f} rho_c={r['rho']:+9.4f}"
              f" sing={r['sing']:.1e} degen3={r['degen3']:.3f}")


if __name__ == "__main__":
    main(float(sys.argv[1]), float(sys.argv[2]))
