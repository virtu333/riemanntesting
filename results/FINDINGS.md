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

## Session 3 (Step 2: statistics to T = 5000)

Machinery. `scan2.py`: both constituents from one batch of the four shared
Hurwitz zetas; truncation N = 0.55 t validated to ~1e-12 on the line and
~2e-11 across sigma in [-1.5, 2.5] up to t = 5000 (5.6x over the 1.1 t
default); vectorized bisection; missed-zero insurance = slope-alternation
check + dip detector + dedup. Cross-validation: reproduces all 397 + 397
constituent zeros at T <= 500 to 9e-12. Result: 5799 + 5801 zeros to 5000.
Winding reuses verified hunt2 with f on the dual evaluator (regression: the
anchor window and a full-N t = 4003 window reproduce exactly).

Pipeline counts for (500, 5000]: 4148 mutual-nearest pairs -> 2940 in-phase
+ 1208 anti-phase -> 821 cancellation events -> 286 retreat events (35%).
Winding 264 merged retreat windows: **237 off-line mirror pairs, zero
bookkeeping mismatches in all 309 windows** (controls included). Controls:
30 non-retreat anti-phase events -> 0 off-line; 49/286 survivors (17%,
Lehmer-like pairs of Z_DH), each re-verified with hardened windows.

Completeness audit (blind wall-to-wall winding of [1000,1100], [2000,2100],
[3000,3100], [4000,4100], [4900,5000] — 11% of the range, no screening):
23 off-line zeros vs 20 in the screened catalog => **screen completeness
87%, 3 misses**, all at t > 3000, all at anti-phase events scored just past
the no-retreat threshold (Theta'_min/omega = +0.16 .. +0.36) with the
LARGEST kappas in the data (0.63-0.72). Failure mode identified at
t = 3059.5: a constituent-cluster configuration (three wave-1 zeros within
0.9 units, mean spacing 0.78, with amplitude collapse |Z'| 17.9 -> 2.2) that
the pairwise two-wave model does not describe (fidelity 0.16 there vs 0.05
typical). The largest-excursion events preferentially come from this cluster
channel — it matters for tails. Practical fix for future runs: widen the
screen to Theta'_min/omega < 0.4.
  (Process note: the audit itself initially had a silent coverage gap whose
  location was CORRELATED with the physics — wrong-sign regions widen the
  line-zero gap that edge-snapping depends on. Fixed by deriving all segment
  boundaries from one whole-block line scan with asserted continuity; the
  gap had hidden exactly one off-line pair, confirmed real by mpmath at 30
  digits.)

Catalog: `data/offline_zeros_5000.json` — 240 zeros in (500, 5000] (3
flagged screen_miss from the audit), plus the 10 known at T <= 500: 250
total. Estimated true count in (500, 5000]: ~272 +- 20 (87% completeness on
an 11% sample).

**Excursion law (the main quantitative result).** The bifurcation picture
(off-line pair = wrong-sign extremum of height h ~ kappa (A+B), curvature
~ (A+B) omega^2, so Im t = sqrt(2h/|Z''|)) predicts with NO free parameters:
    sigma - 1/2 = sqrt(2 kappa) / omega(t),  omega = 0.5 log(5t/2pi).
Measured on 237 events: median actual/predicted = 1.017, r = 0.84; free
two-variable fit gives sigma-1/2 ~ 1.16 kappa^0.533 omega^-0.88 (R^2 = 0.58)
— both exponents land on the predicted (0.5, -1). Log-sd scatter ~44%.

**Growth of N_off(T) — finite-range statement only.** 250 zeros in
(0, 5000], counts per 1250-unit bin = [38, 70, 66, 76]: a constant rate is
REJECTED (chi2 = 13.6/3 dof, p ~ 0.004); over the measured range the rate is
locally better fit by log^2 t (chi2 = 1.8) than log t (chi2 = 4.8). Same
ordering by Poisson logL (Delta = 8.1 and 1.9). Mechanism decomposition per
bin: anti-phase event density scales like omega (ratio 1.16 vs predicted
1.17); retreat fraction flat at 0.34-0.36; hit rate mildly rising 0.76 ->
0.87. IMPORTANT: a log^2 rate CANNOT be asymptotic — the total zero count is
~ (2/pi) T omega(T) ~ T log T, and Bombieri-Hejhal implies N_off = o(T log T)
for such combinations, so the observed super-linear rate is necessarily
transient. The correct reading: over T <= 5000 the off-line event rate grows
faster than constant; the eventual crossover is itself a research question
(the coefficient-family experiment below probes it). Incompleteness skews
high-t, so the finite-range rate is at least as steep as observed.

Caveats: kappa here is the model kappa from the local fit (5% typical
reconstruction error); completeness is measured, not assumed; all "laws" are
empirical fits at T <= 5000.
