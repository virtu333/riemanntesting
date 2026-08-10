# Findings log

## Session 1 (T <= 500 census + mechanism)
- 398 zeros in strip, 378 on line, 10 off-line mirror pairs; no mismatches.
- Necessity law: 230 in-phase collisions -> 0 off-line; all 10 off-line pairs
  from 74 anti-phase collisions.
- corr(excursion, kappa) = 0.75; corr(excursion, gap) = 0.83 (n=10).
- Off-line pair = wrong-sign extremum of Z_DH (failed Lehmer pair).
- Open: the second variable deciding lift-off among anti-phase collisions.

## Session 2 (Step 1b: the second variable)

- Local two-wave model per collision, from constituent data only: PCHIP phase
  through the exact zero phases (slope-sign parity fixes phi mod 2pi), PCHIP
  log-amplitude through |Z'|/phi' at zeros plus |Z| at the extrema between
  them. Median reconstruction error 4.8% of window max (worst 31%). A linear
  phase fit does NOT work: spacing fluctuations are order-1 in phase and its
  intercept error at the collision reaches pi/2.
- Retreat criterion. For W = A e^{i phi1} + B e^{i phi2} (Z_DH = Re W),
  Theta' = [A^2 w1 + B^2 w2 + AB(w1+w2)cos dphi + (A'B - AB')sin dphi]/E^2.
  Constant amplitudes with common carrier give Theta' = w > 0 identically —
  a pure two-wave cancellation can never lose zeros. Zero loss needs phase
  retreat (min Theta' < 0), driven by envelope drift (A'B - AB') and/or
  frequency mismatch weighted by (A - B).
- Event-level result at T <= 500 (74 anti-phase collisions cluster into 56
  cancellation events; adjacent collisions share overlapping windows):
  **13 retreat events = ALL 10 off-line events + 3 survivors.**
  In-phase control: 1/230 collisions retreat (t = 238.9, marginal at -0.5 w).
  Identical confusion matrices for fit windows K = 3, 4, 5.
- The 3 survivors are Lehmer-like pairs of Z_DH (deep cancellation, stayed on
  line, very tight zero pair):
    t = 285.93  pair gap 0.26 msp, extremum +0.208
    t = 457.54  pair gap 0.17 msp, extremum -0.078
    t = 475.55  pair gap 0.27 msp, extremum +0.218
- Roadmap hypothesis FALSIFIED: no slow-variable pair tried separates the
  branch perfectly — not (kappa, nu = dlog(A/B)/dt), not Theta'_min, not even
  the EXACT envelope depth: survivor t=457.54 cancels to 0.028 (A+B units),
  deeper than off-line event t=114.18 at 0.053. Deep cancellation + anti-phase
  + retreat is still not sufficient.
- The decisive object is the sign of a discriminant: off-line zeros of f are
  exactly complex zeros of Z_DH(t) at Im t = -(sigma - 1/2), so lift-off is
  the local model's zero pair going complex. Deciding that sign requires model
  fidelity comparable to the cancellation depth kappa; at ~5% fidelity the
  deep events are undecidable. The analytic-continuation predictor got 6/10
  events, predicted excursions r = 0.57 vs actual (n = 6): right order of
  magnitude, unreliable sign. This sensitivity is intrinsic (bifurcation
  margin), not an implementation defect.
- Necessity chain now: anti-phase (74/304 collisions) -> retreat (13/56
  events) -> off-line (10). Retreat needs line data only, so it is a cheap
  Step 2 screen on top of the anti-phase screen. n = 13 events: every
  quantitative statement above is a hypothesis, not a law.
