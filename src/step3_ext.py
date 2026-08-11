"""Extended-census analysis to T = 20000 (census-only above 5000).

The census below 5000 is validated against wound catalogs at r = 1 and
r = 0.1; above 5000 it stands on the same machinery (spot-validations are
run separately).  Content:
 1. sd(log r*) in eight t-bins against the Selberg conjecture
    var(log r*) = a * log log t + b  (the lever is now ~11% in sd).
 2. Collapse of Psi_T(lambda) with genuinely different sd_T
    (T = 5000 / 10000 / 20000).
 3. T*(r) into the small-r tail (r down to 1e-3).
 4. Normality of log r* at n ~ 2900; duality symmetry at scale.
Writes ../results/phi_ext.json.
"""
import json
import numpy as np
import census_analysis
census_analysis.TAGS = ["c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7"]
from census_analysis import load
from math import erfc, erf, sqrt


def signed_count(t, r, side, T, rv):
    m = (r > 0) & (r <= rv) & (t <= T)
    return int(np.sum(side[m] == 1) - np.sum(side[m] == -1))


def main():
    c = load()
    t, r, side = c["t"], c["r"], c["side"]
    pos = r > 0
    lr, tp = np.log(r[pos]), t[pos]
    out = {}
    print(f"extended census: {len(t)} folds to T={t.max():.0f}, "
          f"r*>0: {pos.sum()}")

    # ---- 1. variance vs log log t ----
    edges = np.array([0, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000])
    rows = []
    print("\nvar(log r*) vs log log t (Selberg: linear):")
    for a, b in zip(edges, edges[1:]):
        m = (tp > a) & (tp <= b)
        if m.sum() < 30:
            continue
        v = float(np.var(lr[m]))
        ll = float(np.log(np.log(np.median(tp[m]))))
        rows.append((ll, v, int(m.sum())))
        print(f"  ({a:5.0f},{b:5.0f}] n={m.sum():4d}  var={v:.3f}  "
              f"loglog(tmed)={ll:.3f}")
    ll_, v_, n_ = map(np.array, zip(*rows))
    wfit = np.polyfit(ll_, v_, 1, w=np.sqrt(n_))
    print(f"  weighted linear fit: var = {wfit[0]:.2f} * loglog t "
          f"{wfit[1]:+.2f}   (Selberg-type slope should be positive O(1))")
    out["var_fit"] = dict(slope=float(wfit[0]), intercept=float(wfit[1]),
                          bins=[dict(loglog=a, var=b, n=int(cn))
                                for a, b, cn in rows])

    # ---- 2. collapse with a longer lever ----
    lg = np.linspace(0.2, 2.5, 12)
    mats, sds = [], {}
    for T in (5000, 10000, 20000):
        sd = float(np.std(lr[tp <= T]))
        sds[T] = sd
        n1 = signed_count(t, r, side, T, 1.0)
        rg = np.exp(np.linspace(np.log(1e-3), 0, 50))
        lam = np.log(1 / rg) / sd
        psi = np.array([signed_count(t, r, side, T, rv) for rv in rg]) / n1
        mats.append(np.interp(lg, lam[::-1], psi[::-1]))
        print(f"\nT={T}: sd={sd:.3f}  N_off(T,1)={n1}")
    mats = np.array(mats)
    spread = mats.max(0) - mats.min(0)
    print("collapse (Psi at common lambda; spread across T=5k/10k/20k):")
    for i, l in enumerate(lg):
        print(f"  lambda={l:4.2f}  Psi={mats.mean(0)[i]:.3f}  "
              f"spread={spread[i]:.3f}")
    out["collapse"] = dict(lambda_grid=lg.tolist(),
                           psi=[m.tolist() for m in mats],
                           sds=sds, spread=spread.tolist())

    # ---- 3. T*(r) tail ----
    births = pos & (side == 1)
    print("\nT*(r) with the extended range:")
    tstars = []
    for rv in (0.1, 0.03, 0.01, 0.003, 0.001):
        m = births & (r <= rv)
        meas = float(t[m].min()) if m.any() else None
        n20 = signed_count(t, r, side, 20000, rv)
        tstars.append(dict(r=rv, T_star=meas, N_off_20000=n20))
        print(f"  r={rv:6.3f}  T*={meas}  N_off(20000, r)={n20}")
    out["Tstar"] = tstars

    # ---- 4. normality + duality at scale ----
    z = (lr - lr.mean()) / lr.std()
    zs = np.sort(z)
    cdf = 0.5 * (1 + np.array([erf(x / sqrt(2)) for x in zs]))
    emp = (np.arange(len(zs)) + 0.5) / len(zs)
    ks = float(np.max(np.abs(cdf - emp)))
    print(f"\nlog r* at scale (n={len(lr)}): mean {lr.mean():+.4f} "
          f"(duality: 0), skew {float(np.mean(z**3)):+.3f}, "
          f"kurt {float(np.mean(z**4)-3):+.3f}, KS {ks:.4f} "
          f"(crit {1.36/np.sqrt(len(lr)):.4f})")
    out["normality"] = dict(n=int(len(lr)), mean=float(lr.mean()),
                            skew=float(np.mean(z**3)),
                            kurt=float(np.mean(z**4) - 3), ks=ks)
    json.dump(out, open("../results/phi_ext.json", "w"), indent=1)
    print("saved ../results/phi_ext.json")


if __name__ == "__main__":
    main()
