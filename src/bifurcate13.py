"""Fold censuses for mod-13 two-wave families (universality test).

Families (indices into mod13.waves order [w2, w10, w4, w8, w6]):
    A : (0, 1)  conjugate pair, order-6 characters   (DH analog)
    B : (2, 3)  conjugate pair, order-3 characters   (DH analog)
    M : (0, 2)  mixed pair, non-conjugate
    C : (4, 0)  real-character vertex + complex-pair wave

Usage: python3 bifurcate13.py FAM T1 T2   (FAM in A,B,M,C;
                                           writes ../data/census13_<FAM>.npz)
"""
import sys
import numpy as np
import bifurcate
import mod13

FAMS = dict(A=(0, 1), B=(2, 3), M=(0, 2), C=(4, 0))

if __name__ == "__main__":
    fam, T1, T2 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    i, j = FAMS[fam]
    bifurcate.PAIR = mod13.pair_factory(i, j)
    bifurcate.OMEGA = lambda t: 0.5 * np.log(13 * t / (2 * np.pi))
    rows = []
    for t1 in np.arange(T1, T2, 12.5):
        t2 = min(t1 + 12.5, T2)
        rows += bifurcate.census_window(t1, t2)
        print(f"[{t1:7.1f},{t2:7.1f}] {fam}: {len(rows)} points", flush=True)
    np.savez(f"../data/census13_{fam}.npz",
             **{k: np.array([r[k] for r in rows]) for k in rows[0].keys()})
    npos = sum(1 for r in rows if 0 < r["r"] <= 1)
    print(f"saved ../data/census13_{fam}.npz: {len(rows)} points, "
          f"{npos} with 0 < r* <= 1")
