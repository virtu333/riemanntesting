"""Off-line zero hunt for Davenport-Heilbronn: full-strip phase winding in log space.

Usage: python3 hunt2.py T1 T2     (windows of 12.5; writes hunt_<T1>_<T2>.json)

Method per window [t1, t2]:
  1. line_count: sign changes of the real Z-function on the critical line (on-line zeros)
  2. count_strip: winding number of completed F around box [-1.5, 2.5] x [t1, t2]
     computed purely from phases (loggamma + angle(f)) so nothing overflows.
     t1 > 0 keeps Gamma poles (s = -1, -3, ...) and trivial zeros (real axis) outside.
  3. discrepancy = total - line = 2 x (mirror pairs off the line, since zeros come in
     (sigma, 1-sigma) pairs); localize the right-half member by t-then-sigma bisection
     of plain-f winding boxes, polish with complex Newton.

Engineering pitfalls learned the hard way:
  - NEVER put a contour edge at sigma = 0.5 + tiny: on-line zeros sit on the boundary
    and blow up adaptive refinement. Use the full symmetric strip instead.
  - Chunk all vectorized evals (~1500 pts): the (N_terms x N_pts) Hurwitz matrix
    otherwise eats GBs at large t.
  - Adaptive refinement: if any per-step |d phase| > 1.5 rad, shrink steps 3x, retry.
"""
import numpy as np
import scipy.special as sp
import json
from dh_core import f_dh

LOG5PI = np.log(5 / np.pi)

def f_chunked(s, chunk=1500):
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    out = np.empty_like(s)
    for i in range(0, len(s), chunk):
        out[i:i + chunk] = f_dh(s[i:i + chunk])
    return out

def phase_F(s):
    """arg of completed F, computed without overflow."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    lg = sp.loggamma((s + 1) / 2)
    return lg.imag + ((s + 1) / 2).imag * LOG5PI + np.angle(f_chunked(s))

def wind_closed(pts, phases):
    d = np.diff(phases)
    d = (d + np.pi) % (2 * np.pi) - np.pi
    return d.sum() / (2 * np.pi), np.max(np.abs(d))

def rect_pts(s1, s2, t1, t2, dv, dh_):
    e1 = complex(s1, t1) + np.arange(0, s2 - s1, dh_)
    e2 = complex(s2, t1) + 1j * np.arange(0, t2 - t1, dv)
    e3 = complex(s2, t2) - np.arange(0, s2 - s1, dh_)
    e4 = complex(s1, t2) - 1j * np.arange(0, t2 - t1, dv)
    return np.concatenate([e1, e2, e3, e4, [complex(s1, t1)]])

def count_strip(t1, t2, s1=-1.5, s2=2.5, dv=0.05, dh_=0.02):
    for attempt in range(4):
        pts = rect_pts(s1, s2, t1, t2, dv, dh_)
        w, mx = wind_closed(pts, phase_F(pts))
        if mx < 1.5 and abs(w - round(w)) < 0.15:
            return int(round(w))
        dv /= 3; dh_ /= 3
    return int(round(w))

def count_box_f(s1, s2, t1, t2, dv=0.04, dh_=0.02):
    """Zeros of plain f in a box away from the critical line (right half)."""
    for attempt in range(4):
        pts = rect_pts(s1, s2, t1, t2, dv, dh_)
        w, mx = wind_closed(pts, np.angle(f_chunked(pts)))
        if mx < 1.5 and abs(w - round(w)) < 0.15:
            return int(round(w))
        dv /= 3; dh_ /= 3
    return int(round(w))

def line_count(t1, t2, dt=0.01):
    ts = np.arange(t1, t2, dt)
    s = 0.5 + 1j * ts
    lg = sp.loggamma((s + 1) / 2)
    ph = lg.imag + (ts / 2) * LOG5PI
    Z = (np.exp(1j * ph) * f_chunked(s)).real
    return int(np.sum(np.signbit(Z[:-1]) != np.signbit(Z[1:])))

def localize_right(t1, t2, n_expect):
    """Zeros of f with sigma in (0.501, 2.5): bisect t, then sigma, Newton polish.
    Caveat: an off-line zero with sigma - 0.5 < 1e-3 would be missed; the
    line-vs-strip bookkeeping mismatch flags that case if it ever occurs."""
    found = []
    def rec(a, b, sl, sr, n):
        if n == 0:
            return
        if (b - a) < 0.02 and (sr - sl) < 0.02:
            z = complex(0.5 * (sl + sr), 0.5 * (a + b))
            for _ in range(40):
                h = 1e-6
                fz = f_chunked([z])[0]
                dfz = (f_chunked([z + h])[0] - f_chunked([z - h])[0]) / (2 * h)
                step = fz / dfz
                z -= step
                if abs(step) < 1e-10:
                    break
            found.append((z.real, z.imag, n))
            return
        if (b - a) >= (sr - sl):
            m = 0.5 * (a + b)
            n1 = count_box_f(sl, sr, a, m)
            rec(a, m, sl, sr, n1); rec(m, b, sl, sr, n - n1)
        else:
            m = 0.5 * (sl + sr)
            n1 = count_box_f(sl, m, a, b)
            rec(a, b, sl, m, n1); rec(a, b, m, sr, n - n1)
    n = count_box_f(0.501, 2.5, t1, t2)
    rec(t1, t2, 0.501, 2.5, n)
    return n, found

if __name__ == "__main__":
    import sys
    T1, T2 = float(sys.argv[1]), float(sys.argv[2])
    windows = np.arange(T1, T2, 12.5)
    results = []
    for t1 in windows:
        t2 = t1 + 12.5
        nl = line_count(t1, t2)
        nt = count_strip(t1, t2)
        noff = nt - nl
        row = {"t1": t1, "t2": t2, "line": nl, "total": nt, "off": noff, "zeros": []}
        if noff > 0:
            nr, zs = localize_right(t1, t2, noff // 2)
            row["zeros"] = [(a, b) for a, b, _ in zs]
            check = "OK" if 2 * nr == noff else f"MISMATCH right={nr}"
            print(f"[{t1:6.1f},{t2:6.1f}] line={nl:3d} total={nt:3d} off={noff} {check} "
                  f"zeros={[(round(a,4),round(b,4)) for a,b in row['zeros']]}", flush=True)
        else:
            print(f"[{t1:6.1f},{t2:6.1f}] line={nl:3d} total={nt:3d} off=0", flush=True)
        results.append(row)
    fn = f"hunt_{int(T1)}_{int(T2)}.json"
    json.dump(results, open(fn, "w"))
    print("saved", fn)
