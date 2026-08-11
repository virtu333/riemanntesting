"""Validate the birth diagram away from r = 1: wind the family at r = 0.1.

Two checks in [500, 2000]:
 1. Every census-predicted arc spanning r = 0.1 (folds with r* <= 0.1,
    clustered within 6 mean spacings) is wound with f_{0.1}; the localized
    off-line pairs must match the signed fold count per window.
 2. A wall-to-wall blind strip [1000, 1100] at r = 0.1 (no census input):
    its total off-line count must equal the census signed count there.

Usage: python3 step3_rslice.py       (writes ../results/rslice_0.1.json)
"""
import json
import numpy as np
import hunt2
from census_analysis import load
from census_verify import set_r, G_r

R = 0.1
T1, T2 = 500.0, 2000.0
STRIP = (1000.0, 1100.0)


def line_zeros_r(a, b, dt=0.01):
    ts = np.arange(a, b, dt)
    Z = G_r(ts, R)
    idx = np.where(np.signbit(Z[:-1]) != np.signbit(Z[1:]))[0]
    return ts[idx] + 0.5 * dt


def snap(edge, zt):
    k = np.searchsorted(zt, edge)
    if k == 0 or k == len(zt):
        return edge
    return 0.5 * (zt[k - 1] + zt[k])


def wind_r(a, b):
    zt = line_zeros_r(a - 3.0, b + 3.0)
    a2, b2 = snap(a, zt), snap(b, zt)
    nline = int(np.sum((zt > a2) & (zt < b2)))
    ntot = hunt2.count_strip(a2, b2)
    off = ntot - nline
    row = dict(a=float(a2), b=float(b2), line=nline, total=ntot, off=off,
               zeros=[], mismatch="")
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
    set_r(R)
    c = load()
    t, r, side = c["t"], c["r"], c["side"]
    sel = (r > 0) & (r <= R) & (t > T1) & (t <= T2)
    folds = sorted(zip(t[sel], side[sel], r[sel]))
    print(f"folds with r* <= {R} in ({T1},{T2}]: {len(folds)} "
          f"(signed count {sum(s for _, s, _ in folds)})")
    # cluster folds within 6 mean spacings
    clusters = []
    for tt, ss, rr in folds:
        msp = np.pi / (0.5 * np.log(5 * tt / (2 * np.pi)))
        if clusters and tt - clusters[-1][-1][0] < 6 * msp:
            clusters[-1].append((tt, ss, rr))
        else:
            clusters.append([(tt, ss, rr)])
    out = dict(r=R, windows=[], strip=None)
    ok = True
    for cl in clusters:
        tmid = 0.5 * (cl[0][0] + cl[-1][0])
        msp = np.pi / (0.5 * np.log(5 * tmid / (2 * np.pi)))
        expect = int(sum(s for _, s, _ in cl))
        row = wind_r(cl[0][0] - 4 * msp, cl[-1][0] + 4 * msp)
        row["expect_pairs"] = expect
        match = (row["off"] == 2 * expect) and not row["mismatch"]
        ok &= match
        out["windows"].append(row)
        print(f"[{row['a']:8.2f},{row['b']:8.2f}] folds={len(cl)} "
              f"expect {expect} pair(s): off={row['off']} "
              f"zeros={[(round(s,4), round(tt,4)) for s, tt in row['zeros']]} "
              f"{'OK' if match else 'MISMATCH ' + row['mismatch']}", flush=True)
    # blind strip
    zt = line_zeros_r(STRIP[0] - 5, STRIP[1] + 5)
    bounds = [snap(x, zt) for x in np.arange(STRIP[0], STRIP[1] + 1e-9, 5.0)]
    tot_off, zs_strip, mm = 0, [], 0
    for a2, b2 in zip(bounds, bounds[1:]):
        nline = int(np.sum((zt > a2) & (zt < b2)))
        ntot = hunt2.count_strip(a2, b2)
        off = ntot - nline
        tot_off += max(off, 0)
        if off > 0:
            nr, zs = hunt2.localize_right(a2, b2, off // 2)
            zs_strip += [(float(s), float(tt)) for s, tt, _ in zs]
            if 2 * nr != off:
                mm += 1
        print(f"strip [{a2:8.2f},{b2:8.2f}] off={off}", flush=True)
    selS = (r > 0) & (r <= R) & (t > STRIP[0]) & (t <= STRIP[1])
    census_strip = int(np.sum(side[selS] == 1) - np.sum(side[selS] == -1))
    out["strip"] = dict(range=STRIP, off_pairs=int(tot_off // 2),
                        census=census_strip, zeros=zs_strip, mismatches=mm)
    print(f"\nstrip {STRIP}: wound {tot_off//2} pair(s), census predicts "
          f"{census_strip}, zeros={zs_strip}, mismatches={mm}")
    print("ALL WINDOWS MATCH" if ok else "SOME WINDOWS MISMATCH")
    json.dump(out, open("../results/rslice_0.1.json", "w"), indent=1)


if __name__ == "__main__":
    main()
