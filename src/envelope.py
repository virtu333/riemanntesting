"""Step 1b: the second variable deciding lift-off at anti-phase collisions.

Local two-wave model with SLOW parameters fitted from constituent zeros only:
    Z_DH(t) ~ A(t) cos phi1(t) + B(t) cos phi2(t) = Re W(t),
    W = A e^{i phi1} + B e^{i phi2},  E = |W|,  Theta = arg W.
Instantaneous frequency of the sum:
    Theta' = [A^2 w1 + B^2 w2 + A B (w1+w2) cos(dphi) + (A'B - AB') sin(dphi)] / E^2.
With constant amplitudes and w1 = w2 this is w > 0 identically: a pure two-wave
sum with common carrier NEVER loses zeros, no matter how deep the cancellation.
Zero loss (a wrong-sign extremum = off-line lift-off) requires phase retreat
Theta' < 0, driven by (a) frequency mismatch weighted by (A - B), and/or
(b) the envelope-drift term (A'B - AB') sin(dphi), which can dominate near deep
anti-phase cancellation where E^2 is small.  Hypothesis under test: the sign of
min_t Theta'(t) over the collision window separates the 10 off-line outcomes
from the 64 on-line ones among the 74 anti-phase collisions (and predicts
uniformly on-line for the 230 in-phase controls).

Fitting (per collision, K zeros each side per wave, constituent data ONLY):
  - phase: zeros of A cos phi sit at phi = pi/2 + m pi with sign(Z') = -(-1)^m;
    consecutive zeros advance phi by pi, slope-sign parity fixes phi mod 2pi.
    Monotone PCHIP interpolation through the exact zero phases (a LINEAR fit is
    far too rigid: spacing fluctuations are order-1 in phase, and its intercept
    error at the collision point can reach pi/2, destroying dphi).
  - amplitude: at a zero, |Z'| = A * phi', so A_k = |Z'(gamma_k)| / phi'(gamma_k);
    PCHIP through log A_k gives A(t), A'(t).

Exact-side diagnostic (NOT used for prediction): fine scan of Z_DH in the
window; sign changes within +-0.75 mean spacings of the deepest cancellation
(on-line collisions keep their 2 zeros there; lifted ones lose both), and the
exact oscillator envelope sqrt(Z^2 + (Z'/w)^2) minimum.

Usage: python3 envelope.py [K]    (needs g1.npy/g2.npy and data/offline json;
                                   writes ../results/envelope_T500.json;
                                   K = zeros each side in the fits, default 4)
"""
import numpy as np, json, os, sys
from scipy.interpolate import PchipInterpolator
from collisions import Z_single, zeros_of

H = 1e-4          # slope stencil
K = 4             # constituent zeros each side used in the slow fits
NGRID = 400       # model-evaluation grid points per window
DT_EXACT = 0.005  # exact Z_DH scan step


def slopes(g, conj):
    return (Z_single(g + H, conj) - Z_single(g - H, conj)) / (2 * H)


def mutual_pairs(g1, g2):
    j = np.searchsorted(g2, g1); j = np.clip(j, 1, len(g2) - 1)
    out = []
    for i, g in enumerate(g1):
        jj = j[i] if abs(g2[j[i]] - g) < abs(g2[j[i] - 1] - g) else j[i] - 1
        k = np.searchsorted(g1, g2[jj]); k = np.clip(k, 1, len(g1) - 1)
        kk = k if abs(g1[k] - g2[jj]) < abs(g1[k - 1] - g2[jj]) else k - 1
        if kk == i:
            out.append((i, jj))
    return out


def fit_wave(g, d, idx, conj):
    """Local slow model of one constituent around zero index idx.
    Returns (phi, logA, tlo, thi): PCHIP interpolants for phase and log-amplitude
    plus their validity range, or None if too few zeros / parity break.
    Amplitude nodes: |Z'|/phi' at the zeros PLUS direct envelope samples |Z| at
    the extrema between them (phi = m pi, solved on the phase interpolant)."""
    lo, hi = max(0, idx - K), min(len(g), idx + K + 1)
    if hi - lo < 4:
        return None
    gz, dz = g[lo:hi], d[lo:hi]
    # slope-sign parity -> absolute phase mod 2pi: phi_k = pi/2 + m_k pi,
    # m_k = k + p with p in {0,1} chosen so sign(dz) = -(-1)^{m_k}
    k = np.arange(hi - lo)
    p = 0 if np.sign(dz[0]) == -1 else 1
    if not np.all(np.sign(dz) == -(-1.0) ** (k + p)):
        return None  # missed zero / parity break in window: bail, don't guess
    phiv = np.pi / 2 + (k + p) * np.pi
    phi = PchipInterpolator(gz, phiv)
    # phase speed at each zero from neighboring spacings (PCHIP endpoint
    # derivatives can clamp to zero, which a division cannot survive)
    wz = np.pi / np.gradient(gz)
    # extrema: phi = m pi between consecutive zeros, solved on the interpolant
    text = np.empty(len(gz) - 1)
    for m in range(len(gz) - 1):
        tgt = phiv[m] + np.pi / 2
        a, b = gz[m], gz[m + 1]
        for _ in range(40):
            c = 0.5 * (a + b)
            if phi(c) < tgt:
                a = c
            else:
                b = c
        text[m] = 0.5 * (a + b)
    aext = np.abs(Z_single(text, conj))
    nodes = np.concatenate([gz, text])
    vals = np.concatenate([np.abs(dz) / wz, aext])
    o = np.argsort(nodes)
    logA = PchipInterpolator(nodes[o], np.log(np.maximum(vals[o], 1e-12)))
    return phi, logA, gz[0], gz[-1]


def local_cubic(wv, tc):
    """Analytic local model of one wave near tc: cubic through the 4 nearest
    exact zero phases; LS cubic through the 8 nearest amplitude nodes.
    Returns poly coeffs (phi, logA) in x = t - tc, plus the node-span radius."""
    phi, logA = wv[0], wv[1]
    gz = phi.x
    near = np.sort(np.argsort(np.abs(gz - tc))[:4])
    x = gz[near] - tc
    cphi = np.polyfit(x, phi(gz[near]), 3)      # interpolation: 4 pts, deg 3
    an = logA.x
    nearA = np.sort(np.argsort(np.abs(an - tc))[:8])
    xa = an[nearA] - tc
    clA = np.polyfit(xa, logA(an[nearA]), 3)
    span = min(-x[0], x[-1], -xa[0], xa[-1])
    return cphi, clA, span


def model_zeros(cub1, cub2, tc, msp, seeds_im=(0.0, 0.1, 0.3)):
    """Complex zeros of M(t) = A1 cos(phi1) + A2 cos(phi2) near tc via damped
    Newton.  M is real on the real axis, so zeros come in conjugate pairs;
    canonical representatives have Im >= 0.  Off-line zeros of f at
    sigma = 1/2 + a are exactly zeros of Z_DH(t) at Im t = -a, so a conjugate
    pair of model zeros predicts lift-off with excursion |Im t|.
    Zeros outside the fit-node span are extrapolation artifacts: rejected."""
    (cp1, ca1, sp1), (cp2, ca2, sp2) = cub1, cub2
    span = min(sp1, sp2)
    dp1, da1 = np.polyder(cp1), np.polyder(ca1)
    dp2, da2 = np.polyder(cp2), np.polyder(ca2)

    def M(x):
        A1, A2 = np.exp(np.polyval(ca1, x)), np.exp(np.polyval(ca2, x))
        return (A1 * np.cos(np.polyval(cp1, x))
                + A2 * np.cos(np.polyval(cp2, x)))

    def dM(x):
        p1, p2 = np.polyval(cp1, x), np.polyval(cp2, x)
        A1, A2 = np.exp(np.polyval(ca1, x)), np.exp(np.polyval(ca2, x))
        return (A1 * (np.polyval(da1, x) * np.cos(p1)
                      - np.polyval(dp1, x) * np.sin(p1))
                + A2 * (np.polyval(da2, x) * np.cos(p2)
                        - np.polyval(dp2, x) * np.sin(p2)))

    zeros = []
    res = np.linspace(-0.9 * msp, 0.9 * msp, 7)
    for xr in res:
        for xi in seeds_im:
            z, converged = complex(xr, xi), False
            for _ in range(60):
                fz, d = M(z), dM(z)
                if not (np.isfinite(fz) and np.isfinite(d)) or d == 0:
                    break
                step = fz / d
                if abs(step) > 0.4 * msp:            # damp: stay local
                    step *= 0.4 * msp / abs(step)
                z -= step
                if abs(z.real) > 2 * msp or abs(z.imag) > 1.0:
                    break                             # wandered out: give up
                if abs(step) < 1e-12:
                    converged = True
                    break
            if converged and abs(z.real) < span and abs(z.imag) < 0.8:
                z = complex(z.real, abs(z.imag))
                if not any(abs(z - u) < 1e-5 for u in zeros):
                    zeros.append(z)
    return sorted(zeros, key=lambda z: abs(z.real))  # offsets from tc


def theta_prime(ts, wv1, wv2):
    """Model Theta'(t), envelope E(t), amplitudes and dphi at times ts."""
    phi1, logA1, _, _ = wv1
    phi2, logA2, _, _ = wv2
    A, B = np.exp(logA1(ts)), np.exp(logA2(ts))
    w1, w2 = phi1.derivative()(ts), phi2.derivative()(ts)
    s1, s2 = logA1.derivative()(ts), logA2.derivative()(ts)
    dphi = phi1(ts) - phi2(ts)
    E2 = A ** 2 + B ** 2 + 2 * A * B * np.cos(dphi)
    num = (A ** 2 * w1 + B ** 2 * w2 + A * B * (w1 + w2) * np.cos(dphi)
           + (s1 - s2) * A * B * np.sin(dphi))
    return num / np.maximum(E2, 1e-30), np.sqrt(np.maximum(E2, 0)), A, B, dphi


def exact_scan(tc, half):
    ts = np.arange(tc - half, tc + half, DT_EXACT)
    Z = Z_single(ts, False) + Z_single(ts, True)
    return ts, Z


def main():
    if os.path.exists("g1.npy"):
        g1, g2 = np.load("g1.npy"), np.load("g2.npy")
    else:
        g1 = zeros_of(False); g2 = zeros_of(True)
        np.save("g1.npy", g1); np.save("g2.npy", g2)
    d1, d2 = slopes(g1, False), slopes(g2, True)
    cat = json.load(open("../data/offline_zeros_T500.json"))
    off_t = np.array([z["t"] for z in cat])
    off_s = np.array([z["sigma"] for z in cat])

    pairs = mutual_pairs(g1, g2)
    mids = np.array([0.5 * (g1[i] + g2[j]) for i, j in pairs])
    # outcome: nearest-collision assignment (dedup, no window double count)
    outcome = np.zeros(len(pairs), dtype=int)
    for t in off_t:
        outcome[np.argmin(np.abs(mids - t))] = 1

    rows, skipped = [], 0
    for n, (i, j) in enumerate(pairs):
        tc = mids[n]
        wv1 = fit_wave(g1, d1, i, False)
        wv2 = fit_wave(g2, d2, j, True)
        if wv1 is None or wv2 is None:
            skipped += 1
            rows.append(dict(mid=float(tc), sgn=int(np.sign(d1[i] * d2[j])),
                             out=int(outcome[n]), fit="failed"))
            continue
        wbar = 0.5 * float(wv1[0].derivative()(tc) + wv2[0].derivative()(tc))
        msp = np.pi / wbar
        tlo = max(wv1[2], wv2[2], tc - msp)
        thi = min(wv1[3], wv2[3], tc + msp)
        ts = np.linspace(tlo, thi, NGRID)
        thp, E, A, B, dphi = theta_prime(ts, wv1, wv2)
        imin = int(np.argmin(thp))
        icanc = int(np.argmin(E / (A + B)))
        Ac, Bc = A[icanc], B[icanc]
        # model fidelity vs exact Z_DH at 7 points across the window
        tchk = np.linspace(tlo, thi, 7)
        zex = Z_single(tchk, False) + Z_single(tchk, True)
        _, _, Am, Bm, _ = theta_prime(tchk, wv1, wv2)
        zmod = (Am * np.cos(wv1[0](tchk)) + Bm * np.cos(wv2[0](tchk)))
        fid = float(np.max(np.abs(zmod - zex)) / np.max(np.abs(zex)))
        s1c = float(wv1[1].derivative()(ts[icanc]))
        s2c = float(wv2[1].derivative()(ts[icanc]))
        row = dict(
            mid=float(tc), gap=float(abs(g1[i] - g2[j])),
            sgn=int(np.sign(d1[i] * d2[j])), out=int(outcome[n]),
            sig=float(off_s[np.argmin(np.abs(off_t - tc))] - 0.5)
                if outcome[n] else None,
            dw=float(wv1[0].derivative()(ts[icanc]) - wv2[0].derivative()(ts[icanc])),
            nu=float(s1c - s2c),                       # d/dt log(A/B) at cancellation
            rho=float(min(Ac, Bc) / max(Ac, Bc)),
            kappa_model=float(E[icanc] / (Ac + Bc)),
            dphi_c=float(dphi[icanc] % (2 * np.pi)),
            thmin=float(thp[imin]), thmin_norm=float(thp[imin] / wbar),
            t_thmin=float(ts[imin]), t_canc=float(ts[icanc]),
            fid=fid,
        )
        # discriminant predictor: complex zeros of the local analytic model
        tcanc = row["t_canc"]
        zs = model_zeros(local_cubic(wv1, tcanc), local_cubic(wv2, tcanc),
                         tcanc, msp)
        loc = [z for z in zs if abs(z.real) < 0.7 * msp]
        cplx = [z for z in loc if z.imag > 1e-3]
        row["pred_off"] = int(len(cplx) > 0)
        row["pred_exc"] = float(max(z.imag for z in cplx)) if cplx else 0.0
        row["pred_t"] = float(tcanc + cplx[0].real) if cplx else None
        row["n_real_loc"] = len(loc) - len(cplx)
        # exact diagnostic: sign changes near the cancellation point
        if row["sgn"] < 0:
            tse, Z = exact_scan(row["t_canc"], 1.5 * msp)
            near = np.abs(tse - row["t_canc"]) < 0.75 * msp
            zn = Z[near]
            row["nzero_exact"] = int(np.sum(np.signbit(zn[:-1]) != np.signbit(zn[1:])))
            dZ = np.gradient(Z, DT_EXACT)
            osc = np.sqrt(Z ** 2 + (dZ / wbar) ** 2)
            row["Emin_exact"] = float(np.min(osc[near]) / (Ac + Bc))
        rows.append(row)

    os.makedirs("../results", exist_ok=True)
    suffix = "" if K == 4 else f"_K{K}"
    json.dump(rows, open(f"../results/envelope_T500{suffix}.json", "w"), indent=1)

    ok = [r for r in rows if r.get("fit") != "failed"]
    anti = [r for r in ok if r["sgn"] < 0]
    inph = [r for r in ok if r["sgn"] > 0]
    print(f"pairs={len(pairs)} fitted={len(ok)} skipped={skipped} "
          f"anti={len(anti)} in={len(inph)}")
    for name, grp in (("anti", anti), ("in", inph)):
        for o in (0, 1):
            sub = [r for r in grp if r["out"] == o]
            neg = sum(1 for r in sub if r["thmin"] < 0)
            po = sum(1 for r in sub if r["pred_off"])
            print(f"  {name}-phase out={o}: n={len(sub)}  thmin<0: {neg}  "
                  f"pred_off: {po}")
    offr = [r for r in ok if r["out"]]
    print("\npredicted vs actual excursion (the 10 off-line events):")
    for r in offr:
        print(f"  mid={r['mid']:8.3f} actual={r['sig']:.4f} "
              f"pred={r['pred_exc']:.4f} pred_t={r['pred_t']}")
    print("median model fidelity (max |Z_model-Z|/max|Z| per window): "
          f"{np.median([r['fid'] for r in ok]):.3f}  "
          f"worst: {max(r['fid'] for r in ok):.3f}")
    print("\nanti-phase collisions sorted by thmin_norm "
          "(OFF = produced an off-line pair):")
    for r in sorted(anti, key=lambda r: r["thmin_norm"]):
        tag = "OFF" if r["out"] else "   "
        nz = r.get("nzero_exact", "?")
        print(f"  {tag} mid={r['mid']:8.3f} thmin/w={r['thmin_norm']:+8.4f} "
              f"kap={r['kappa_model']:.3f} nu={r['nu']:+.4f} dw={r['dw']:+.4f} "
              f"rho={r['rho']:.3f} nz={nz} "
              f"Emin={r.get('Emin_exact', float('nan')):.3f} fid={r['fid']:.2f}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        K = int(sys.argv[1])
    main()
