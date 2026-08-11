"""Step 2, stage D: assemble the off-line catalog and test the scalings.

Inputs: step2_events.json (screen), step2_wind.json (winding), the T<=500
catalog.  Outputs: data/offline_zeros_5000.json (new zeros in (500, 5000]),
event outcomes (off-line vs survivor) for every retreat event, and the
statistics the roadmap asks for:
  - excursion sigma-1/2 vs cancellation depth kappa (correlation + power fit),
  - predicted excursion (Step 1b discriminant) vs actual, out-of-sample,
  - N_off(T) growth: compare c*T against c*T*log^2(T) by Poisson likelihood
    on equal-count bins,
  - survivor fraction of retreat events vs height.
Controls are checked and reported: every mismatch or control violation prints.
"""
import json
import numpy as np


def main():
    ev = json.load(open("../results/step2_events.json"))
    wd = json.load(open("../results/step2_wind.json"))
    cat0 = json.load(open("../data/offline_zeros_T500.json"))

    # ---- controls and bookkeeping ----
    for name in ("retreat", "ctrl_anti", "ctrl_bg"):
        bad = [r for r in wd[name] if r["mismatch"]]
        noff = sum(r["off"] for r in wd[name])
        print(f"{name:9s}: {len(wd[name])} windows, sum(off)={noff}, "
              f"mismatches={len(bad)}")
        for r in bad:
            print(f"   MISMATCH [{r['a']:.2f},{r['b']:.2f}]: {r['mismatch']}")

    # ---- attribute localized zeros to retreat events ----
    revs = [e for e in ev["events"] if e["retreat"]]
    tcs = np.array([e["t_canc"] for e in revs])
    zeros = [z for r in wd["retreat"] for z in r["zeros"]]
    for e in revs:
        e["off_zeros"] = []
    orphans = []
    for sig, t in zeros:
        k = int(np.argmin(np.abs(tcs - t)))
        if abs(tcs[k] - t) < 3.5 * revs[k]["msp"]:
            revs[k]["off_zeros"].append((sig, t))
        else:
            orphans.append((sig, t))
    if orphans:
        print(f"ORPHAN zeros (no retreat event within 3.5 msp): {orphans}")
    n_off_ev = sum(1 for e in revs if e["off_zeros"])
    multi = [e for e in revs if len(e["off_zeros"]) > 1]
    print(f"\nretreat events: {len(revs)}  -> off-line: {n_off_ev}  "
          f"survivors: {len(revs) - n_off_ev}  "
          f"(events with >1 zero: {len(multi)})")

    newcat = sorted([dict(sigma=s, t=t) for s, t in zeros],
                    key=lambda z: z["t"])
    json.dump(newcat, open("../data/offline_zeros_5000.json", "w"), indent=1)
    print(f"catalog: {len(newcat)} off-line zeros in (500, 5000] "
          f"(+{len(cat0)} known at T<=500)")

    # ---- excursion vs kappa, predicted vs actual ----
    ex, kap, pred = [], [], []
    for e in revs:
        if e["off_zeros"]:
            sig = max(s for s, _ in e["off_zeros"]) - 0.5
            ex.append(sig); kap.append(e["kappa_model"]); pred.append(e["pred_exc"])
    ex, kap, pred = map(np.array, (ex, kap, pred))
    ok = (ex > 0) & (kap > 0)
    r_lin = np.corrcoef(kap[ok], ex[ok])[0, 1]
    r_log = np.corrcoef(np.log(kap[ok]), np.log(ex[ok]))[0, 1]
    a, b = np.polyfit(np.log(kap[ok]), np.log(ex[ok]), 1)
    print(f"\nexcursion vs kappa (n={ok.sum()}): r={r_lin:.3f}  "
          f"log-log r={r_log:.3f}  power fit sigma-1/2 ~ {np.exp(b):.3f} "
          f"* kappa^{a:.3f}")
    both = pred > 0
    if both.sum() > 2:
        r_pred = np.corrcoef(pred[both], ex[both])[0, 1]
        print(f"discriminant pred_exc vs actual (n={both.sum()}): "
              f"r={r_pred:.3f}  (pred>0 rate among off events: "
              f"{both.sum()}/{len(pred)})")
    p0 = np.array([e["pred_exc"] for e in revs if not e["off_zeros"]])
    print(f"pred_exc>0 among survivors: {np.sum(p0 > 0)}/{len(p0)}")

    # ---- N_off(T) growth ----
    ts = np.sort(np.array([z["t"] for z in cat0] +
                          [z["t"] for z in newcat]))
    print(f"\nN_off growth, {len(ts)} zeros total in (0, 5000]:")
    for name, dens in (("T (uniform)", lambda t: np.ones_like(t)),
                       ("T log^2 T", lambda t: np.log(t) ** 2),
                       ("T log T", lambda t: np.log(t))):
        # Poisson log-likelihood of the observed points under rate ~ dens
        tt = np.linspace(1, 5000, 20000)
        norm = np.trapezoid(dens(tt), tt)
        ll = np.sum(np.log(dens(ts) / norm * len(ts)))
        print(f"  rate ~ {name:12s}: logL = {ll:9.3f}")
    half = np.searchsorted(ts, 2500.0)
    print(f"  split check: {half} zeros below t=2500, {len(ts)-half} above")

    json.dump([dict(t_canc=e["t_canc"], kappa_model=e["kappa_model"],
                    thmin_norm=e["thmin_norm"], pred_exc=e["pred_exc"],
                    off_zeros=e["off_zeros"]) for e in revs],
              open("../results/step2_event_outcomes.json", "w"), indent=1)
    print("saved ../results/step2_event_outcomes.json and "
          "../data/offline_zeros_5000.json")


if __name__ == "__main__":
    main()
