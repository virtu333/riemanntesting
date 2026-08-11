"""Analysis of the Wronskian bifurcation census.

Checks, in order of importance:
 1. r = 1 reconstruction (the census must reproduce the DH catalog):
    per t-window, line count changes by -2 at a side=+1 fold and +2 at a
    side=-1 fold as r crosses r*, so
        N_off(W, r) = #(side=+1, 0 < r* <= r) - #(side=-1, 0 < r* <= r).
    Compare to the 250 verified off-line zeros at r = 1 (10 + 240),
    globally and per 500-unit bin.
 2. Fold-law excursion check at r = 1: each catalog zero matched to its
    nearest census fold (0 < r* < 1, complex side toward 1); compare
    C sqrt(1 - r*) with the measured sigma - 1/2 (a LOCAL law used far from
    the fold: agreement should improve as r* -> 1).
 3. r <-> 1/r duality: the chi <-> chibar exchange maps the family to its
    mirror with r* -> 1/r*, so the log r* distribution should be symmetric
    about 0 (statistical check).
 4. Amplitude-matching: in the two-wave picture r* = A/B at the fold;
    correlation of log r* with log(A/B).
 5. The birth diagram: N_off(T, r) on a log-r grid (the primary artifact),
    plus first-escape height T*(r), plus degeneracy-diagnostic tails.
"""
import json
import numpy as np

TAGS = ["c0", "c1", "c2", "c3"]


def load():
    cols = None
    for tag in TAGS:
        d = np.load(f"../data/census_{tag}.npz")
        if cols is None:
            cols = {k: [d[k]] for k in d.files}
        else:
            for k in d.files:
                cols[k].append(d[k])
    cols = {k: np.concatenate(v) for k, v in cols.items()}
    o = np.argsort(cols["t"])
    return {k: v[o] for k, v in cols.items()}


def main():
    c = load()
    t, r, side = c["t"], c["r"], c["side"]
    pos = r > 0
    print(f"census points: {len(t)}  (r*>0: {pos.sum()}, "
          f"r*<0 mirror family: {(~pos).sum()})")
    print(f"r* consistency: median rerr {np.median(c['rerr']):.2e}, "
          f"worst {np.max(c['rerr']):.2e}")

    # ---- 1. r = 1 reconstruction ----
    cat = json.load(open("../data/offline_zeros_T500.json")) \
        + json.load(open("../data/offline_zeros_5000.json"))
    ct = np.array([z["t"] for z in cat])
    cs = np.array([z["sigma"] for z in cat])
    print(f"\ncatalog: {len(cat)} verified off-line zeros at r=1")
    sel = pos & (r <= 1)
    n_pred = int(np.sum(side[sel] == 1) - np.sum(side[sel] == -1))
    print(f"census signed count at r=1, all t: {n_pred}")
    print("per 500-unit bin (census predicted vs catalog):")
    for a in np.arange(0, 5000, 500):
        m = sel & (t > a) & (t <= a + 500)
        npred = int(np.sum(side[m] == 1) - np.sum(side[m] == -1))
        ncat = int(np.sum((ct > a) & (ct <= a + 500)))
        flag = "" if npred == ncat else "   <-- MISMATCH"
        print(f"  ({a:4.0f},{a+500:4.0f}]: census {npred:3d}  "
              f"catalog {ncat:3d}{flag}")

    # ---- 2. fold-law excursion at r=1 ----
    birth = np.where(pos & (r < 1) & (side == 1))[0]
    preds, acts, rstars = [], [], []
    for tt, ss in zip(ct, cs):
        k = birth[np.argmin(np.abs(t[birth] - tt))]
        if abs(t[k] - tt) > 2.0:
            continue
        preds.append(c["C"][k] * np.sqrt(1 - r[k]))
        acts.append(ss - 0.5)
        rstars.append(r[k])
    preds, acts, rstars = map(np.array, (preds, acts, rstars))
    print(f"\nfold-law check at r=1 (n={len(preds)} matched):")
    print(f"  corr(pred, actual) = {np.corrcoef(preds, acts)[0,1]:.3f}, "
          f"median actual/pred = {np.median(acts/preds):.3f}")
    for lo, hi in ((0.6, 1.0), (0.3, 0.6), (0.0, 0.3)):
        m = (rstars >= lo) & (rstars < hi)
        if m.sum() > 3:
            print(f"  r* in [{lo},{hi}): n={m.sum():3d}  "
                  f"median actual/pred = {np.median(acts[m]/preds[m]):.3f}")

    # ---- 3. duality symmetry of log r* ----
    lr = np.log(r[pos])
    print(f"\nlog r* (n={len(lr)}): mean {np.mean(lr):+.4f} "
          f"(duality predicts 0), sd {np.std(lr):.3f}, "
          f"skew {float(((lr-lr.mean())**3).mean()/lr.std()**3):+.3f}")

    # ---- 4. amplitude matching ----
    ab = np.log(c["A"][pos] / c["B"][pos])
    print(f"corr(log r*, log(A/B)) = {np.corrcoef(lr, ab)[0,1]:.3f}  "
          f"(two-wave picture predicts ~1 with slope 1: "
          f"fit slope {np.polyfit(ab, lr, 1)[0]:.3f})")

    # ---- 5. birth diagram + T*(r) + degeneracy ----
    grid = np.exp(np.linspace(np.log(1e-3), np.log(1.0), 25))
    diagram = []
    for rv in grid:
        m = pos & (r <= rv)
        n = int(np.sum(side[m] == 1) - np.sum(side[m] == -1))
        births = t[pos & (r <= rv) & (side == 1)]
        diagram.append(dict(r=float(rv), N_off=n,
                            T_star=float(births.min()) if len(births) else None))
    json.dump(diagram, open("../results/birth_diagram.json", "w"), indent=1)
    print("\nbirth diagram (r, N_off(5000, r), T*(r)):")
    for d in diagram[::4]:
        print(f"  r={d['r']:.4f}  N_off={d['N_off']:4d}  T*={d['T_star']}")
    dg = c["degen"][pos]
    print(f"\ndegeneracy diagnostic |G''|/((A+B) w^2) over r*>0 folds: "
          f"median {np.median(dg):.3f}, 5th pct {np.percentile(dg,5):.4f}, "
          f"1st pct {np.percentile(dg,1):.4f}")
    nsmall = int(np.sum(dg < 0.05))
    print(f"folds with degen < 0.05 (cluster-mediated candidates): {nsmall}")
    json.dump(dict(n_points=int(len(t)), n_pos=int(pos.sum()),
                   reconstruction=n_pred, catalog=len(cat)),
              open("../results/census_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
