"""Step 2, stage B: collision screen over [500, 5000] from line data only.

Loads the scan shards (plus the validated T500 scan for boundary continuity),
builds mutual-nearest cross pairs, classifies in/anti-phase, runs the Step 1b
local two-wave fit + retreat criterion on every anti-phase collision, and
clusters collisions into cancellation events.  No strip winding here: the
output is the candidate list for stage C.

Usage: python3 step2_screen.py       (reads ../data/scan_*.npz,
                                      writes ../results/step2_events.json)
"""
import json
import numpy as np
import envelope
import scan2
from events import cluster_events

# extremum sampling inside fit_wave goes through the fast dual evaluator
envelope.Z_single = lambda ts, conj, chunk=800: scan2.Z_pair(ts)[1 if conj else 0]

TMIN, TMAX = 500.0, 5000.0
TAGS = ["T500", "s1", "s2", "s3", "s4"]


def load_scans():
    g1, d1, g2, d2 = [], [], [], []
    for tag in TAGS:
        s = np.load(f"../data/scan_{tag}.npz")
        g1.append(s["g1"]); d1.append(s["d1"])
        g2.append(s["g2"]); d2.append(s["d2"])
    g1, d1 = np.concatenate(g1), np.concatenate(d1)
    g2, d2 = np.concatenate(g2), np.concatenate(d2)
    o1, o2 = np.argsort(g1), np.argsort(g2)
    return g1[o1], d1[o1], g2[o2], d2[o2]


def main():
    g1, d1, g2, d2 = load_scans()
    print(f"zeros: wave1 {len(g1)}, wave2 {len(g2)} up to T={g1[-1]:.1f}")
    pairs = envelope.mutual_pairs(g1, g2)
    mids = np.array([0.5 * (g1[i] + g2[j]) for i, j in pairs])
    keep = (mids > TMIN) & (mids <= TMAX)
    print(f"mutual-nearest pairs in ({TMIN},{TMAX}]: {int(keep.sum())}")

    rows, n_in, n_skip = [], 0, 0
    for n, (i, j) in enumerate(pairs):
        if not keep[n]:
            continue
        sgn = int(np.sign(d1[i] * d2[j]))
        if sgn > 0:
            n_in += 1
            continue
        tc = mids[n]
        wv1 = envelope.fit_wave(g1, d1, i, False)
        wv2 = envelope.fit_wave(g2, d2, j, True)
        if wv1 is None or wv2 is None:
            n_skip += 1
            continue
        wbar = 0.5 * float(wv1[0].derivative()(tc) + wv2[0].derivative()(tc))
        msp = np.pi / wbar
        tlo = max(wv1[2], wv2[2], tc - msp)
        thi = min(wv1[3], wv2[3], tc + msp)
        ts = np.linspace(tlo, thi, envelope.NGRID)
        thp, E, A, B, dphi = envelope.theta_prime(ts, wv1, wv2)
        imin = int(np.argmin(thp))
        icanc = int(np.argmin(E / (A + B)))
        zs = envelope.model_zeros(
            envelope.local_cubic(wv1, ts[icanc]),
            envelope.local_cubic(wv2, ts[icanc]), ts[icanc], msp)
        cplx = [z for z in zs if abs(z.real) < 0.7 * msp and z.imag > 1e-3]
        rows.append(dict(
            mid=float(tc), gap=float(abs(g1[i] - g2[j])), sgn=sgn,
            nu=float(wv1[1].derivative()(ts[icanc])
                     - wv2[1].derivative()(ts[icanc])),
            dw=float(wv1[0].derivative()(ts[icanc])
                     - wv2[0].derivative()(ts[icanc])),
            kappa_model=float(E[icanc] / (A[icanc] + B[icanc])),
            rho=float(min(A[icanc], B[icanc]) / max(A[icanc], B[icanc])),
            thmin=float(thp[imin]), thmin_norm=float(thp[imin] / wbar),
            t_thmin=float(ts[imin]), t_canc=float(ts[icanc]), msp=float(msp),
            pred_exc=float(max(z.imag for z in cplx)) if cplx else 0.0,
        ))
    print(f"anti-phase fitted: {len(rows)}  in-phase: {n_in}  "
          f"fit failures: {n_skip}")
    evs = cluster_events(rows)
    out = []
    for e in evs:
        r = min(e, key=lambda r: r["thmin"])
        out.append(dict(t_canc=r["t_canc"], msp=r["msp"], n_collisions=len(e),
                        thmin_norm=r["thmin_norm"],
                        retreat=int(r["thmin"] < 0),
                        kappa_model=r["kappa_model"], nu=r["nu"], dw=r["dw"],
                        rho=r["rho"], pred_exc=max(x["pred_exc"] for x in e)))
    n_re = sum(e["retreat"] for e in out)
    print(f"anti-phase events: {len(out)}   retreat events: {n_re}")
    json.dump(dict(summary=dict(n_pairs=int(keep.sum()), n_in=n_in,
                                n_anti=len(rows), n_skip=n_skip,
                                n_events=len(out), n_retreat=n_re),
                   collisions=rows, events=out),
              open("../results/step2_events.json", "w"), indent=1)
    print("saved ../results/step2_events.json")


if __name__ == "__main__":
    main()
