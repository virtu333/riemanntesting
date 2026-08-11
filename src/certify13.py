"""Independent certification of cusp candidates by mpmath Newton (30 dps).

For each candidate (t, r, rho) from the census, run Newton on the full
system (G, G', G'') = 0 in (t, r, rho) with all derivatives computed by
mpmath's high-precision differentiation of the Hurwitz-zeta representation
— fully independent of the census's finite-difference machinery.  Record
convergence, the displacement from census values (= empirical census
error), the final residual, G''' (must be bounded away from zero), and the
Jacobian condition number at the solution.

Usage: python3 certify13.py SHARD NSHARDS  (stratified sample hardcoded;
       writes ../results/cusp_cert_shard<SHARD>.jsonl)
"""
import json
import sys
import numpy as np
import mpmath as mp

mp.mp.dps = 30

IDX = {}
v = 1
for k in range(12):
    IDX[v] = k
    v = (v * 2) % 13


def chi(j, a):
    return mp.e ** (2j * mp.pi * j * IDX[a] / 12)


def L(j, s):
    return 13 ** (-s) * mp.fsum([chi(j, a) * mp.zeta(s, mp.mpf(a) / 13)
                                 for a in range(1, 13)])


def delta(j):
    tau = mp.fsum([chi(j, a) * mp.e ** (2j * mp.pi * a / 13)
                   for a in range(1, 13)])
    return 0.5 * mp.arg(tau / mp.sqrt(13))


D2, D10, D6 = delta(2), delta(10), delta(6)


def wave(j, d, t):
    s = mp.mpc(0.5, t)
    th = mp.im(mp.loggamma(s / 2)) + t / 2 * mp.log(13 / mp.pi)
    return mp.re(mp.e ** (1j * (th - d)) * L(j, s))


def G(t, r, rho):
    return wave(2, D2, t) + r * wave(10, D10, t) + rho * wave(6, D6, t)


def Gp(t, r, rho):
    return mp.diff(lambda x: G(x, r, rho), t)


def Gpp(t, r, rho):
    return mp.diff(lambda x: G(x, r, rho), t, 2)


def certify(t0, r0, rho0, tol=mp.mpf(10) ** -18, iters=8):
    t0, r0, rho0 = mp.mpf(t0), mp.mpf(r0), mp.mpf(rho0)
    tin, rin, rhoin = t0, r0, rho0
    for _ in range(iters):
        F = mp.matrix([G(t0, r0, rho0), Gp(t0, r0, rho0),
                       Gpp(t0, r0, rho0)])
        J = mp.matrix(3, 3)
        for col in range(3):
            for row, fn in enumerate((G, Gp, Gpp)):
                args = [t0, r0, rho0]

                def f1(x, args=args, col=col, fn=fn):
                    a = list(args)
                    a[col] = x
                    return fn(*a)
                J[row, col] = mp.diff(f1, args[col])
        try:
            st = mp.lu_solve(J, F)
        except ZeroDivisionError:
            return dict(converged=False, reason="singular Jacobian")
        t0, r0, rho0 = t0 - st[0], r0 - st[1], rho0 - st[2]
        res = max(abs(G(t0, r0, rho0)), abs(Gp(t0, r0, rho0)),
                  abs(Gpp(t0, r0, rho0)))
        if res < tol:
            break
    g3 = mp.diff(lambda x: G(x, r0, rho0), t0, 3)
    # condition number of the final Jacobian
    Jn = np.array([[float(J[i, k]) for k in range(3)] for i in range(3)])
    s = np.linalg.svd(Jn, compute_uv=False)
    return dict(
        converged=bool(res < tol),
        t=float(t0), r=float(r0), rho=float(rho0),
        dt=float(abs(t0 - tin)), dr=float(abs(r0 - rin)),
        drho=float(abs(rho0 - rhoin)),
        residual=float(res), g3=float(g3),
        cond=float(s[0] / s[-1]))


def stratified_sample(n=30, seed=11):
    c = np.load("../data/cusps13.npz")
    ok = (c["boundary"] == 0) & (c["t"] > 30)
    idx = np.where(ok)[0]
    rng = np.random.default_rng(seed)
    # stratify by t-quintile and |r| magnitude (small/large)
    tq = np.quantile(c["t"][idx], [0.2, 0.4, 0.6, 0.8])
    strata = {}
    for i in idx:
        b = int(np.searchsorted(tq, c["t"][i]))
        m = int(abs(c["r"][i]) > 2)
        strata.setdefault((b, m), []).append(i)
    picks = []
    per = max(1, n // len(strata))
    for k, members in sorted(strata.items()):
        picks += list(rng.choice(members, min(per, len(members)),
                                 replace=False))
    # always include the two headline cusps
    for tt in (1510.764119, 102.967699):
        j = int(np.argmin(np.abs(c["t"] - tt)))
        if j not in picks:
            picks.append(j)
    return [(float(c["t"][i]), float(c["r"][i]), float(c["rho"][i]),
             float(c["g3"][i])) for i in picks]


if __name__ == "__main__":
    shard, nsh = int(sys.argv[1]), int(sys.argv[2])
    sample = stratified_sample()
    fh = open(f"../results/cusp_cert_shard{shard}.jsonl", "w")
    for k, (t0, r0, rho0, g3fd) in enumerate(sample):
        if k % nsh != shard:
            continue
        out = certify(t0, r0, rho0)
        out.update(census_t=t0, census_r=r0, census_rho=rho0,
                   census_g3=g3fd)
        fh.write(json.dumps(out) + "\n")
        fh.flush()
        if out.get("converged"):
            print(f"t={t0:9.3f}: CERT dt={out['dt']:.2e} dr={out['dr']:.2e} "
                  f"drho={out['drho']:.2e} resid={out['residual']:.1e} "
                  f"G'''={out['g3']:+.3f} (fd {g3fd:+.3f}) "
                  f"cond={out['cond']:.1e}", flush=True)
        else:
            print(f"t={t0:9.3f}: FAILED {out.get('reason','')}", flush=True)
    fh.close()
    print("shard done")
