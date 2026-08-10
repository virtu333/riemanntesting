"""Event-level summary of the Step 1b envelope analysis.

Adjacent mutual-nearest collisions whose evaluation windows overlap detect the
same cancellation point; cluster collisions into EVENTS by t_canc (within 0.7)
and evaluate the retreat criterion (min model Theta' < 0) per event.

For retreat events that stayed on the line ("survivors"), scan exact Z_DH and
record the surviving zero pair and the extremum height between them: these are
Lehmer-like pairs of Z_DH, the near-critical objects of the collision picture.

Usage: python3 events.py          (reads ../results/envelope_T500.json,
                                   writes ../results/events_anti_T500.json)
"""
import json
import numpy as np
from collisions import Z_single, omega


def cluster_events(anti):
    evs = []
    for r in sorted(anti, key=lambda r: r["t_canc"]):
        if evs and abs(r["t_canc"] - evs[-1][-1]["t_canc"]) < 0.7:
            evs[-1].append(r)
        else:
            evs.append([r])
    return evs


def survivor_pair(tc):
    """Exact zeros of Z_DH bracketing tc and the extremum between them."""
    ts = np.arange(tc - 2.5, tc + 2.5, 0.002)
    Z = Z_single(ts, False) + Z_single(ts, True)
    sc = np.where(np.signbit(Z[:-1]) != np.signbit(Z[1:]))[0]
    zt = ts[sc]
    i = np.argsort(np.abs(zt - tc))[:2]
    z1, z2 = sorted(zt[i])
    seg = (ts > z1) & (ts < z2)
    h = float(Z[seg][np.argmax(np.abs(Z[seg]))]) if seg.any() else float("nan")
    return float(z1), float(z2), h


def main():
    rows = [r for r in json.load(open("../results/envelope_T500.json"))
            if r.get("fit") != "failed"]
    anti = [r for r in rows if r["sgn"] < 0]
    evs = cluster_events(anti)
    out = []
    for e in evs:
        r = min(e, key=lambda r: r["thmin"])
        ev = dict(t_canc=r["t_canc"], n_collisions=len(e),
                  thmin_norm=r["thmin_norm"], retreat=int(r["thmin"] < 0),
                  off=int(any(x["out"] for x in e)),
                  sig=next((x["sig"] for x in e if x["out"]), None),
                  kappa_model=r["kappa_model"],
                  Emin_exact=r.get("Emin_exact"), nu=r["nu"], dw=r["dw"],
                  pred_exc=max(x["pred_exc"] for x in e))
        if ev["retreat"] and not ev["off"]:
            z1, z2, h = survivor_pair(ev["t_canc"])
            msp = np.pi / omega(ev["t_canc"])
            ev.update(pair=(z1, z2), pair_gap_msp=float((z2 - z1) / msp),
                      extremum=h)
        out.append(ev)
    json.dump(out, open("../results/events_anti_T500.json", "w"), indent=1)

    n_re = sum(e["retreat"] for e in out)
    n_off = sum(e["off"] for e in out)
    caught = sum(e["retreat"] for e in out if e["off"])
    print(f"anti-phase events: {len(out)}  retreat: {n_re}  off-line: {n_off}"
          f"  off caught by retreat: {caught}/{n_off}")
    print("survivors (retreat but on-line):")
    for e in out:
        if e["retreat"] and not e["off"]:
            print(f"  t={e['t_canc']:8.3f} thmin/w={e['thmin_norm']:+7.2f} "
                  f"kap={e['kappa_model']:.3f} pair_gap={e['pair_gap_msp']:.3f}"
                  f" extremum={e['extremum']:+.4f}")


if __name__ == "__main__":
    main()
