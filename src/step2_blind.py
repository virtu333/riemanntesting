"""Step 2 completeness audit: blind contiguous winding, no screening.

Winds every ~5-unit segment of a contiguous block wall-to-wall (edges snapped
to midpoints between Z_DH line zeros; adjacent segments snap to the same
midpoint, so coverage is exact with no gaps or double counting) and localizes
every off-line zero.  Comparing against the screened catalog measures the
retreat screen's miss rate directly — prompted by a background control window
that caught an off-line pair (sigma-1/2 = 0.368, t = 3059.52) at an anti-phase
event the screen scored as no-retreat (a constituent-cluster configuration the
pairwise two-wave model does not describe; model fidelity 0.16 there).

Coverage discipline (learned the hard way — the first version had a coverage
GAP near t = 3059.5 that silently skipped an off-line pair, because per-segment
edge snapping with a small line-scan margin fell back to an unsnapped edge
precisely where a wrong-sign region widens the line-zero gap): ALL segment
boundaries come from ONE whole-block line scan, each boundary is the midpoint
of the zero gap containing its nominal position, consecutive segments share
boundaries exactly, and continuity is asserted.

Usage: python3 step2_blind.py T1 T2 TAG   (writes
           ../results/blind_<TAG>.jsonl, one row per segment)
"""
import json
import sys
import numpy as np
from step2_wind import wind_core, line_zeros, snap  # fast dual-L patch too


def main(t1, t2, tag):
    zt = line_zeros(t1 - 5.0, t2 + 5.0)
    bounds = [snap(x, zt) for x in np.arange(t1, t2 + 1e-9, 5.0)]
    assert all(b2 > b1 for b1, b2 in zip(bounds, bounds[1:])), "bad bounds"
    path = f"../results/blind_{tag}.jsonl"
    fh = open(path, "w")
    prev_b = None
    for a2, b2 in zip(bounds, bounds[1:]):
        assert prev_b is None or a2 == prev_b, "coverage gap!"
        prev_b = b2
        row = wind_core(a2, b2, zt)
        fh.write(json.dumps(row) + "\n")
        fh.flush()
        flag = f"  {row['mismatch']}" if row["mismatch"] else ""
        print(f"[{row['a']:9.3f},{row['b']:9.3f}] line={row['line']:3d} "
              f"total={row['total']:3d} off={row['off']}"
              f" zeros={[(round(s, 4), round(t, 4)) for s, t in row['zeros']]}"
              f"{flag}", flush=True)
    fh.close()
    print(f"coverage [{bounds[0]:.3f}, {bounds[-1]:.3f}] wall-to-wall; saved",
          path)


if __name__ == "__main__":
    main(float(sys.argv[1]), float(sys.argv[2]), sys.argv[3])
