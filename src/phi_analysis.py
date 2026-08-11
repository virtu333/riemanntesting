"""Step 3a: the crossover function Phi of the singular transition.

All inputs are the Wronskian census (exact fold set).  Content:
 1. Psi_T(r) = N_off(T, r) / N_off(T, 1) for several T-prefixes, and the
    collapse test in lambda = log(1/r) / sd_T(log r*): if the amplitude
    log-ratio at folds is the driving Gaussian variable, curves at
    different T should collapse onto one shape, erfc-like in lambda.
 2. Normality of log r* (moments + KS against a fitted Gaussian), overall
    and in t-bins; binned sd against the sqrt(log log t) conjecture
    (declared weakly testable: ~8% predicted variation over this range).
 3. First-escape height T*(r): measured versus a Poisson first-arrival
    model with intensity nu(t) * P_Gauss(log r* <= log r) where nu(t) is
    the measured birth-fold density (per unit t, proportional to omega
    empirically) and the Gaussian uses the measured (mu, sd) of log r*.
 4. The in-phase / anti-phase fold asymmetry (mirror family folds
    outnumber ours ~3:1) and the phase-offset explanation: dphi is
    centered at -2 delta ~ -0.554, so level 0 is closer than level +-pi
    for the fluctuating phase difference; a Gaussian crossing-rate
    estimate with the measured phase sd is compared to the ratio.

Writes ../results/phi_analysis.json; prints the summary.
"""
import json
import numpy as np
from census_analysis import load
from dh_core import DELTA


def signed_count(t, r, side, T, rv):
    m = (r > 0) & (r <= rv) & (t <= T)
    return int(np.sum(side[m] == 1) - np.sum(side[m] == -1))


def main():
    c = load()
    t, r, side = c["t"], c["r"], c["side"]
    pos = r > 0
    lr = np.log(r[pos])
    tp = t[pos]

    out = {}

    # ---- 1. collapse test ----
    Ts = [1250, 2500, 3750, 5000]
    rgrid = np.exp(np.linspace(np.log(3e-3), 0.0, 40))
    curves = {}
    for T in Ts:
        sd = float(np.std(lr[tp <= T]))
        n1 = signed_count(t, r, side, T, 1.0)
        lam = np.log(1 / rgrid) / sd
        psi = np.array([signed_count(t, r, side, T, rv) for rv in rgrid]) / n1
        curves[T] = dict(sd=sd, n1=n1, lam=lam.tolist(), psi=psi.tolist())
        print(f"T={T:5d}: N_off(T,1)={n1:4d}  sd(log r*, t<=T)={sd:.3f}")
    # collapse quality: interpolate each curve onto common lambda grid
    lg = np.linspace(0.2, 2.5, 12)
    mat = []
    for T in Ts:
        lam = np.array(curves[T]["lam"])[::-1]
        psi = np.array(curves[T]["psi"])[::-1]
        mat.append(np.interp(lg, lam, psi))
    mat = np.array(mat)
    spread = np.max(mat, 0) - np.min(mat, 0)
    mean = np.mean(mat, 0)
    print("\ncollapse in lambda (Psi mean, max spread across T):")
    for i, l in enumerate(lg):
        print(f"  lambda={l:4.2f}  Psi={mean[i]:.3f}  spread={spread[i]:.3f}")
    out["collapse"] = dict(lambda_grid=lg.tolist(), psi_mean=mean.tolist(),
                           spread=spread.tolist(), curves=curves)
    # erfc shape: compare Psi_mean with c * Phi_gauss-tail?  The natural
    # normalized-tail model: Psi(lambda) ~ P(Z > lambda - lambda0)/P(Z > -l0)
    from math import erfc
    resid = []
    for l0 in np.linspace(-1.0, 1.0, 81):
        pred = np.array([erfc((l - l0) / np.sqrt(2)) for l in lg]) \
             / erfc(-l0 / np.sqrt(2))
        resid.append((np.mean((pred - mean) ** 2), l0))
    best = min(resid)
    l0 = best[1]
    pred = np.array([erfc((l - l0) / np.sqrt(2)) for l in lg]) \
        / erfc(-l0 / np.sqrt(2))
    print(f"\nbest erfc fit: lambda0={l0:+.3f}, rms resid "
          f"{np.sqrt(best[0]):.4f}")
    for i, l in enumerate(lg):
        print(f"  lambda={l:4.2f}  Psi={mean[i]:.3f}  erfc-model={pred[i]:.3f}")
    out["erfc_fit"] = dict(lambda0=float(l0), rms=float(np.sqrt(best[0])))

    # ---- 2. normality of log r* ----
    z = (lr - lr.mean()) / lr.std()
    skew = float(np.mean(z ** 3))
    kurt = float(np.mean(z ** 4) - 3)
    zs = np.sort(z)
    from math import erf
    cdf = 0.5 * (1 + np.array([erf(x / np.sqrt(2)) for x in zs]))
    emp = (np.arange(len(zs)) + 0.5) / len(zs)
    ks = float(np.max(np.abs(cdf - emp)))
    print(f"\nlog r* normality (n={len(lr)}): skew {skew:+.3f}, "
          f"excess kurtosis {kurt:+.3f}, KS {ks:.4f} "
          f"(KS 95% crit ~ {1.36/np.sqrt(len(lr)):.4f})")
    print("binned sd(log r*) vs sqrt(log log t):")
    bins = [(0, 1250), (1250, 2500), (2500, 3750), (3750, 5000)]
    binsd = []
    for a, b in bins:
        m = (tp > a) & (tp <= b)
        s = float(np.std(lr[m]))
        ll = np.sqrt(np.log(np.log(0.5 * (a + b) + 1)))
        binsd.append(dict(bin=[a, b], n=int(m.sum()), sd=s, sqrtloglog=ll))
        print(f"  ({a:4d},{b:4d}] n={m.sum():3d}  sd={s:.3f}  "
              f"sqrt(loglog tmid)={ll:.3f}  ratio={s/ll:.3f}")
    out["normality"] = dict(skew=skew, kurt=kurt, ks=ks, binned=binsd)

    # ---- 3. T*(r): measured vs Poisson first-arrival ----
    births = pos & (side == 1)
    tb, lrb = t[births], np.log(r[births])
    mu, sdb = float(lrb.mean()), float(lrb.std())
    # birth density per unit t ~ prop to omega(t): calibrate the constant
    w = lambda x: 0.5 * np.log(5 * x / (2 * np.pi))
    tt = np.linspace(10, 5000, 2000)
    cw = len(tb) / np.trapezoid(w(tt), tt)   # births per unit (omega-weighted)
    from math import sqrt
    print(f"\nT*(r): measured first birth vs Poisson-Gaussian model "
          f"(mu={mu:.3f}, sd={sdb:.3f}):")
    tstar_rows = []
    for rv in (0.3, 0.1, 0.03, 0.01, 0.003):
        m = births & (r <= rv)
        meas = float(t[m].min()) if m.any() else None
        ptail = 0.5 * np.array([erfc(-(np.log(rv) - mu) / (sdb * sqrt(2)))])[0]
        lam = cw * ptail * np.cumsum(w(tt)) * (tt[1] - tt[0])
        k = np.searchsorted(lam, np.log(2))
        pred = float(tt[k]) if k < len(tt) else None
        tstar_rows.append(dict(r=rv, measured=meas, predicted_median=pred))
        print(f"  r={rv:6.3f}  measured T*={meas}  model median T*={pred}")
    out["Tstar"] = tstar_rows

    # ---- 4. in/anti asymmetry ----
    n_anti, n_in = int(pos.sum()), int((~pos).sum())
    # phase-difference sd inferred from the level populations of a Gaussian
    # stationary process crossing levels at distance d0 = 2 delta and pi-2delta
    d0, dpi = 2 * DELTA, np.pi - 2 * DELTA
    ratio = n_in / n_anti
    # crossing rate ~ exp(-d^2 / 2 s^2): solve for s from the ratio
    s2 = (dpi ** 2 - d0 ** 2) / (2 * np.log(ratio))
    print(f"\nfold asymmetry: in-phase {n_in} vs anti-phase {n_anti} "
          f"(ratio {ratio:.2f}); Gaussian crossing model with offset "
          f"2 delta = {d0:.3f} implies phase sd ~ {np.sqrt(s2):.3f} rad")
    out["asymmetry"] = dict(n_in=n_in, n_anti=n_anti, ratio=ratio,
                            implied_phase_sd=float(np.sqrt(s2)))

    json.dump(out, open("../results/phi_analysis.json", "w"), indent=1)
    print("\nsaved ../results/phi_analysis.json")


if __name__ == "__main__":
    main()
