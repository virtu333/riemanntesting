"""Step 2, stage C: strip winding at retreat events + out-of-sample controls.

Windows: every retreat event +-3 mean spacings (overlaps merged), plus two
control families that test the necessity chain out-of-sample:
  - NCTRL_ANTI random NON-retreat anti-phase events (prediction: 0 off-line),
  - NCTRL_BG random 5-unit background windows (prediction: 0 off-line unless
    they happen to contain a retreat event; those are excluded).
Window edges are snapped to midpoints between consecutive Z_DH line zeros so
no on-line zero sits on the contour.  Per window: line count (Z1 + Z2 on the
fast dual evaluator; the identity Z_DH = Z_chi + Z_chibar is exact), full-strip
winding (verified hunt2 machinery, conservative N), off = total - line, and
localization + Newton polish of right-half zeros when off > 0.  Bookkeeping
(total = line + 2 x localized) is checked per window and mismatches reported,
never smoothed over.

Usage: python3 step2_wind.py         (reads ../results/step2_events.json,
                                      writes ../results/step2_wind.json)
"""
import json
import numpy as np
import scan2
import hunt2

RNG = np.random.default_rng(42)
NCTRL_ANTI = 30
NCTRL_BG = 15
HALF_MSP = 3.0


def line_zeros(a, b, dt=0.01):
    ts = np.arange(a, b, dt)
    Z1, Z2 = scan2.Z_pair(ts)
    Z = Z1 + Z2
    idx = np.where(np.signbit(Z[:-1]) != np.signbit(Z[1:]))[0]
    return ts[idx] + 0.5 * dt


def snap(edge, zt):
    """Midpoint of the consecutive-zero gap containing edge."""
    k = np.searchsorted(zt, edge)
    if k == 0 or k == len(zt):
        return edge
    return 0.5 * (zt[k - 1] + zt[k])


def wind_window(a, b):
    zt = line_zeros(a - 1.0, b + 1.0)
    a2, b2 = snap(a, zt), snap(b, zt)
    nline = int(np.sum((zt > a2) & (zt < b2)))
    ntot = hunt2.count_strip(a2, b2)
    off = ntot - nline
    row = dict(a=a2, b=b2, line=nline, total=ntot, off=off, zeros=[],
               mismatch="")
    if off < 0 or off % 2:
        row["mismatch"] = f"odd/negative off={off}"
        return row
    if off > 0:
        nr, zs = hunt2.localize_right(a2, b2, off // 2)
        row["zeros"] = [(float(s), float(t)) for s, t, _ in zs]
        if 2 * nr != off:
            row["mismatch"] = f"right={nr} vs off={off}"
    return row


def main():
    data = json.load(open("../results/step2_events.json"))
    evs = data["events"]
    retreat = [e for e in evs if e["retreat"]]
    quiet = [e for e in evs if not e["retreat"]]
    print(f"{len(retreat)} retreat events, {len(quiet)} quiet anti events")

    # merge overlapping retreat windows
    ivs = sorted([e["t_canc"] - HALF_MSP * e["msp"],
                  e["t_canc"] + HALF_MSP * e["msp"]] for e in retreat)
    merged = []
    for iv in ivs:
        if merged and iv[0] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], iv[1])
        else:
            merged.append(iv)

    ctrl_anti = RNG.choice(len(quiet), min(NCTRL_ANTI, len(quiet)),
                           replace=False)
    ctrl_anti_iv = [[quiet[k]["t_canc"] - HALF_MSP * quiet[k]["msp"],
                     quiet[k]["t_canc"] + HALF_MSP * quiet[k]["msp"]]
                    for k in ctrl_anti]
    bg = []
    while len(bg) < NCTRL_BG:
        a = float(RNG.uniform(500, 4995))
        iv = [a, a + 5.0]
        if not any(iv[0] < m[1] and m[0] < iv[1] for m in merged):
            bg.append(iv)

    out = dict(retreat=[], ctrl_anti=[], ctrl_bg=[])
    for name, group in (("retreat", merged), ("ctrl_anti", ctrl_anti_iv),
                        ("ctrl_bg", bg)):
        for a, b in group:
            row = wind_window(a, b)
            out[name].append(row)
            flag = f"  {row['mismatch']}" if row["mismatch"] else ""
            print(f"{name:9s} [{a:7.1f},{b:7.1f}] line={row['line']:3d} "
                  f"total={row['total']:3d} off={row['off']}"
                  f" zeros={[(round(s, 4), round(t, 4)) for s, t in row['zeros']]}"
                  f"{flag}", flush=True)
    json.dump(out, open("../results/step2_wind.json", "w"), indent=1)
    print("saved ../results/step2_wind.json")


if __name__ == "__main__":
    main()
