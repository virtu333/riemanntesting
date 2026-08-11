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

Speed: hunt2's f evaluation is monkeypatched to the dual-L path (both L's from
the four shared Hurwitz zetas, truncation N = 0.55 t validated to ~2e-11 over
sigma in [-1.5, 2.5] up to t = 5000) — ~5x per point over the conservative
default; work is sharded over cores.  The anchor window and a full-N spot
window are regression targets (see step2_wind_check.py history in FINDINGS).

Usage: python3 step2_wind.py SHARD NSHARDS   (writes incremental
           ../results/step2_wind_shard<SHARD>.jsonl, one JSON row per window)
       python3 step2_wind.py merge           (merges shards ->
           ../results/step2_wind.json)
"""
import json
import sys
import numpy as np
import scan2
import hunt2

# fast f: e^{-i delta} L(s,chi) + e^{i delta} L(s,chibar), one dual evaluation
from dh_core import DELTA


def f_fast(s):
    L1, L2 = scan2.L_pair(s)
    return np.exp(-1j * DELTA) * L1 + np.exp(1j * DELTA) * L2


hunt2.f_dh = f_fast

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


def build_windows():
    """Deterministic full window list: (kind, a, b) for retreat + controls."""
    data = json.load(open("../results/step2_events.json"))
    evs = data["events"]
    retreat = [e for e in evs if e["retreat"]]
    quiet = [e for e in evs if not e["retreat"]]

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
    return ([("retreat", a, b) for a, b in merged]
            + [("ctrl_anti", a, b) for a, b in ctrl_anti_iv]
            + [("ctrl_bg", a, b) for a, b in bg])


def main(shard, nshards):
    wins = build_windows()
    path = f"../results/step2_wind_shard{shard}.jsonl"
    fh = open(path, "w")
    for k, (kind, a, b) in enumerate(wins):
        if k % nshards != shard:
            continue
        row = wind_window(a, b)
        row["kind"] = kind
        fh.write(json.dumps(row) + "\n")
        fh.flush()
        flag = f"  {row['mismatch']}" if row["mismatch"] else ""
        print(f"{kind:9s} [{a:7.1f},{b:7.1f}] line={row['line']:3d} "
              f"total={row['total']:3d} off={row['off']}"
              f" zeros={[(round(s, 4), round(t, 4)) for s, t in row['zeros']]}"
              f"{flag}", flush=True)
    fh.close()
    print(f"shard {shard}/{nshards} done -> {path}")


def merge(nshards):
    out = dict(retreat=[], ctrl_anti=[], ctrl_bg=[])
    for i in range(nshards):
        for line in open(f"../results/step2_wind_shard{i}.jsonl"):
            row = json.loads(line)
            out[row.pop("kind")].append(row)
    for v in out.values():
        v.sort(key=lambda r: r["a"])
    json.dump(out, open("../results/step2_wind.json", "w"), indent=1)
    print("merged", {k: len(v) for k, v in out.items()},
          "-> ../results/step2_wind.json")


if __name__ == "__main__":
    if sys.argv[1] == "merge":
        merge(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    else:
        main(int(sys.argv[1]), int(sys.argv[2]))
