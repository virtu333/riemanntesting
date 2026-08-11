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

## Session 4 (coefficient family and the Wronskian bifurcation census)

Family: G_r(t) = Z1(t) + r Z2(t), r real; DH is r = 1; r = 0 is the genuine
Euler-product vertex L(s, chi) (on-line zeros there = GRH, not a theorem —
none detected in the verified range).  docs/family_lemma.md proves the
Schwarz self-dual FE and the (sigma, t) <-> (1-sigma, t) mirror symmetry for
every real r (verified numerically to ~1e-12, src/family.py).

**Census (`src/bifurcate.py`).** Double zeros G_r = G_r' = 0 are the zeros
of the r-independent Wronskian H = Z1 Z2' - Z1' Z2, with critical
coefficient r* = -Z1/Z2 = -Z1'/Z2' (median consistency 2.4e-6).  In the
two-wave picture H ~ A B omega sin(dphi): the Step-1 anti-phase necessity
law is the sign structure of H (dphi = pi gives r* > 0, our branch; dphi = 0
gives r* < 0, the mirror family Z1 - |r| Z2).  Census over (0, 5000]:
2816 folds, 678 with r* > 0.

**Exact reconstruction of the catalog.** Line count changes by -/+ 2 as r
crosses a fold, so N_off(W, r) = #(side=+1, r* <= r) - #(side=-1, r* <= r)
with no arc-pairing needed.  The census signed count at r = 1 reproduced the
verified catalog EXACTLY in every 500-unit bin after (a) discovering three
new off-line zeros at census-vs-catalog divergence sites — (0.8978,
2442.5303), (0.5236, 2923.4527), (0.8535, 4392.3764), all screen misses,
one the smallest excursion on record — and (b) accounting one boundary-flux
case (the zero at t = 4999.57 whose birth fold lies just above T = 5000).
Catalog now **253 verified zeros in (0, 5000]**
(`data/offline_zeros_5000.json`).  The census supersedes the retreat screen
as the discovery instrument: it is exact (no envelope model), complete up to
H-scan resolution, and cheaper than winding.  True screen completeness in
hindsight: 237/243 = 97.5% (the audit's 87% was small-sample fluctuation).

**Fold normal form verified.**  |sigma - 1/2| = |2 Z2(t*)/G''(t*)|^{1/2}
|r - r*|^{1/2}.  At delta = |r - r*| = 0.01: measured/predicted = 0.997,
0.994, 0.985 (complex side) and real-side gap/(2 C sqrt(delta)) = 1.003,
1.034, 1.036, at folds t* = 85.9, 2452.5, 2923.6.  Drift at delta = 0.1
(ratios 0.94-0.82) is the expected higher-order correction.  At r = 1 the
extrapolated fold law improves as r* -> 1 (median actual/pred 0.897 for
r* in [0.6, 1)) — it is a local law, used far from the fold.

**Structure of the fold set.**
- Duality: chi <-> chibar maps r* -> 1/r*; measured mean log r* = +0.007,
  skew +0.017 (predicted 0).  sd(log r*) = 1.756 (sqrt(log log 5000) = 1.46;
  same scale, ~20% larger).
- Amplitude matching is EXACT at folds: corr(log r*, log(A/B)) = 1.000,
  slope 1.000 (oscillator envelopes at t*): r* = A/B.  The epsilon-scaling
  question reduces to the statistics of log(A/B) at Wronskian zeros.
- Prefactor bridge: two-wave reduction predicts C ~ sqrt(2/(1+r*))/omega;
  measured median C/prediction = 0.825, corr 0.925.  At r* -> 1 this is
  1/omega and recovers the empirical sigma - 1/2 = sqrt(2 kappa)/omega law
  (kappa(r=1) = |1 - r*|/(1 + r*) in the same reduction).
- **No cusp candidates in this one-parameter family**: degeneracy
  diagnostic |G''|/((A+B) omega^2) has 1st percentile 0.13, no fold below
  0.05.  The "cluster-mediated" events are ordinary folds of the exact
  function; the earlier pairwise-model failures were failures of the
  envelope approximation, not a second bifurcation class.  A genuine cusp
  (G = G' = G'' = 0) needs two coefficient parameters — that is the mod-13
  simplex's job.

**Birth diagram (`results/birth_diagram.json`).**  N_off(5000, r) =
252, 180, 67, 8, 2, 0 at r = 1, 0.32, 0.1, 0.032, 0.01, 0.0032; first
escape height T*(r) = 85.9, 177, 241, 922, 2563, none.  This is the
finite-height singular transition: the off-line population switches off
through the tail of the log(A/B) distribution as r -> 0.  Next: predict
N_off(r) from the measured log(A/B) distribution at Wronskian zeros
(erfc-shape test), r-slice blind validation, then the mod-13 two-parameter
simplex for universality and the cusp search.

## Session 5 (the crossover function Phi; `src/phi_analysis.py`)

- Psi_T(r) = N_off(T, r)/N_off(T, 1) for T = 1250/2500/3750/5000 plotted
  against lambda = log(1/r)/sd_T(log r*) collapses onto one curve (max
  spread 0.09 for lambda < 0.7, < 0.02 for lambda > 1.2, transient 0.17 at
  lambda ~ 0.9 from the N = 38 lowest-T curve).  CAVEAT stated up front:
  sd_T varies only 1.65 -> 1.76 over these heights, so the collapse mostly
  reflects a single fixed shape; the T-normalization is weakly tested at
  these heights.  The SHAPE is the content.
- Shape: a one-parameter erfc model fits to rms 0.068 (lambda0 = +0.60)
  but the measured tail decays FASTER than erfc (Psi = 0.005 vs model
  0.040 at lambda = 2.5).  Consistently, log r* is Gaussian by KS
  (0.039 < 0.052 crit at n = 678, skew +0.02) but platykurtic (excess
  kurtosis -0.47): light tails.  Both observations match sub-Gaussian
  large deviations of log|L| at finite height — the CLT (Selberg) regime
  holds in the bulk, decays faster in the tail.  A quantitative
  large-deviation correction is the natural theory target.
- Binned sd(log r*): 1.67/1.64/1.88/1.79 across four t-bins — no resolved
  growth (sqrt(log log t) predicts +7% across the range; noise is ~ +-6%).
  The T = 20000 census extension is running to lengthen this lever.
- T*(r): measured (177/241/922/2563 at r = 0.3/0.1/0.03/0.01) vs a crude
  Poisson first-arrival model with Gaussian tail (52/112/509/4508): order
  of magnitude only.  Model lacks death-fold correction and clustering;
  recorded as a baseline, not a fit.
- Fold asymmetry: in-phase (mirror family) folds outnumber anti-phase
  3.15 : 1.  A Gaussian level-crossing model with the phase offset
  2 delta = 0.554 (dphi is centered at -2 delta; level 0 is closer than
  +-pi) implies phase-difference sd ~ 1.67 rad — the same scale as
  sd(log r*) = 1.76, as Selberg statistics predict for log-modulus vs
  argument fluctuations.  A coherent two-Gaussian picture of the fold set.
- **r-slice validation at r = 0.1 (`src/step3_rslice.py`)**: all 19
  census-predicted arcs in (500, 2000] wound to exactly one off-line pair
  each (19/19, zero bookkeeping mismatches; excursions 0.52-0.72 in
  sigma), and a wall-to-wall blind strip [1000, 1100] with no census input
  found exactly the 2 pairs the census predicts there.  The birth diagram
  is verified away from r = 1.

### Extension to T = 20000 (`src/step3_ext.py`; census-only above 5000)

13404 folds to T = 20000, 3512 in the r > 0 branch.  IMPORTANT framing:
above T = 5000 all N_off values are CENSUS-IMPLIED, not wound-verified —
only three high-t folds were spot-wound (t* = 6885, 12839, 18823; each
produced exactly the predicted off-line pair, positions within 0.04 in t,
excursion ratios 0.82/0.98/0.81 at |1 - r*| = 0.04-0.11).  The census
methodology is validated below 5000 (exact reconstruction) and spot-checked
above; "1267 verified pairs at T = 20000" would overstate it.

- **Variance growth consistent with the Selberg scale**: var(log r*) grows
  monotonically 2.72 -> 3.56 across eight t-bins.  With slope FIXED at the
  Selberg-difference prediction (var = log log t + c; log r* = log A -
  log B, each ~ (1/2) log log t), bins 2-8 fit with c ~ 1.2, residuals
  <= 0.19; only the lowest bin (low-t transient) deviates.  The free-fit
  slope (2.66) is not meaningful on this lever length; the claim is
  "monotone growth consistent with the predicted scale", not a
  demonstrated normalization — the Selberg layer remains a conditional
  asymptotic interpretation.
- **Collapse with a real lever**: sd_T = 1.756 / 1.816 / 1.859 and
  N_off(T, 1) = 252 / 567 / 1267 at T = 5000 / 10000 / 20000; Psi_T(lambda)
  spread <= 0.032 over the full lambda range (was 0.17 with the T = 1250
  curve).  The crossover shape is stable in T.
- **Three decades of transition**: N_off(20000, r) = 1267 / 390 / 81 / 18 /
  5 / 2 at r = 1 / 0.1 / 0.03 / 0.01 / 0.003 / 0.001; first-escape heights
  T*(0.003) = 5000.07, T*(0.001) = 15705.1.
- log r* at n = 3512: mean -0.0008 (duality), skew +0.045, excess kurtosis
  -0.172, KS 0.0235 vs crit 0.0229 — Gaussian bulk with mild light-tail
  deviations, borderline at this n; the platykurtosis shrank with height
  (was -0.47 at T <= 5000), consistent with CLT convergence.
- Growth at r = 1 with the validated census: N_off doubles as 252 -> 567 ->
  1267 for T doublings; the T-ratios (2.25, 2.23) sit between T log T
  (2.16, 2.15) and T log^2 T (2.34, 2.31) — the log-power remains
  unresolved, now with 5x the data and no completeness caveat.

## Session 6 (mod-13 five-wave simplex: universality and the cusp hunt)

Machinery (`src/mod13.py`, validated to 1e-12/5e-13 like mod 5, all five
Euler vertices clean): the even primitive characters mod 13 — real chi_6
plus conjugate pairs (chi_2, chi_10) and (chi_4, chi_8) — share one gamma
factor, giving FIVE real rotated waves spanning a single Schwarz-self-dual
simplex.  Fold censuses at T <= 2500 for four two-wave families
(`src/bifurcate13.py`, injectable-pair refactor of the census; mod-5
anchor regression passes):

  family                     folds  r*>0  mean log r*   sd    r*=A/B  degen p1
  A = (chi2, chi10) conj      1364   234    -0.089    1.717   1.0000   0.181
  B = (chi4, chi8)  conj      1553   551    +0.033    1.864   1.0000   0.080
  M = (chi2, chi4)  mixed     1476   129    +0.196    1.594   0.9998   0.221
  C = (chi6, chi2)  real+cx   1430   198    -0.051    1.612   1.0000   0.119
  mod-5 baseline (t<=2500)    1276   291    -0.018    1.649   1.0000   0.148

Universality findings (with statistical corrections, session 7):
- r* = A/B holds to corr 1.0000 in every family — but this is an
  ALGEBRAIC IDENTITY for our oscillator-envelope definition (at a fold
  Z1 = -r* Z2 and Z1' = -r* Z2', so A = |r*| B exactly), NOT independent
  evidence.  Correctly stated: a lemma verified numerically, whose content
  is that log r* IS the log-amplitude-ratio entering the Selberg
  statistics, by construction.
- Duality symmetry of log r* — RESOLVED with two additional mixed
  families (M2 = (chi2, chi8), n = 428, mean +0.023; M3 = (chi10, chi4),
  n = 458, mean -0.003): the pooled mixed sample (n = 1015) has mean
  +0.033, bootstrap CI [-0.07, +0.14], sign-flip permutation p = 0.53.
  ALL families, conjugate and mixed, are consistent with a symmetric
  log r* distribution; the M family's +0.196 was a fluctuation.  The
  correct statement: symmetry does not require the chi <-> chibar
  exchange argument — the amplitude statistics of any two same-conductor
  L-functions are exchangeable, and the data reflect exactly that.
  (Two earlier phrasings of this point — "the exception proves the
  mechanism", then "consistent with absence of constraint" — are both
  superseded.)
- No near-degenerate folds in ANY one-parameter family (degen 1st
  percentile 0.08-0.22): cusps require two parameters, universally.
- The anti-phase fold fraction (8.7%-35% across families) tracks the
  root-number phase offset via the level-crossing picture: family B
  (offset 2 delta_4 = 1.32, nearly equidistant between levels 0 and pi)
  has the largest anti-phase share; M (effective offset 0.39) the least.
- Two-wave fold-prefactor reduction C ~ sqrt(2/(1+r*))/omega: median
  ratio 0.63-1.18, log-corr 0.88-0.93 across families — same quality as
  mod 5.

### The cusp census and the second universality class (`src/cusp13.py`,
`src/cusp13_verify.py`)

For the three-wave family G = Z1 + r Z2 + rho Z3 (waves chi_2, chi_10,
chi_6), a genuine cusp G = G' = G'' = 0 is a zero of the 3x3 Wronskian
D(t) = det[(Z_i, Z_i', Z_i'')] — a pole-free scalar; (1, r_c, rho_c) is
the null vector.  Census over (0, 2500]: **1467 cusp CANDIDATES** (SVD
residuals ~1e-10; only one in the all-positive coefficient cone).  Cusps
are abundant with two coefficients, absent with one — the codimension
count made flesh.  Candidate status: roots of the numerically
differentiated Wronskian refined with the same machinery; promotion to
"certified cusps" requires independent high-precision derivatives,
step-size/stencil sensitivity, rank-2 checks, and G''' bounds on a
stratified sample (queued).  Coefficient-space framing: the five waves
span a real projective family; the POSITIVE cone containing the
Euler-product vertices holds exactly one detected cusp candidate
(t_c = 1510.764, r_c = 3.475, rho_c = 2.637) — the showcased cusp
(r_c = 6.27, rho_c = -1.94) is a valid Schwarz-self-dual member outside
that cone.  The positive-cone cusp is queued for direct verification and,
if clean, becomes the headline example.

Verification at the cusp t_c = 102.9677 (r_c = 6.2716, rho_c = -1.9398):
- Multiplicity bookkeeping: at the cusp the triple zero counts 3 in the
  strip winding and 1 as a line sign change — off = 2 without any
  off-line zero (engineering note: multiplicity, not lift-off).
- **Generic rays give the CUBE-ROOT law**: excursion exponents 0.310
  (rho-ray) and 0.330 (r-ray) vs predicted 1/3, prefactor
  (sqrt(3)/2)|6 c_0 / G'''|^{1/3} with c_0 the constant unfolding term
  (= drho Z3 or dr Z2), accurate to ~20-30%.  ANY transversal line
  through a cusp sees the cube root — the naive expectation of 1/2 along
  the r-ray is wrong because dr couples to the constant term dr Z2(t_c).
- **The tangent ray (drho = -(Z2/Z3) dr, killing the constant term) gives
  the exact FOLD law**: off = 0 on one side (three real zeros), a complex
  pair on the other with fitted exponent 0.499 vs 1/2 and prefactor
  |6 gamma_eff d / G'''|^{1/2} (gamma_eff = W(Z2,Z3)/Z3) matching to
  <1% over two decades (0.0205/0.0205 ... 0.1179/0.1185).

### Session 7 certification: the positive-cone cusp is the headline example

Both headline cusps independently certified by 30-digit mpmath Newton on
the full system (G, G', G'') = 0 (derivatives via mpmath differentiation
of the Hurwitz representation — independent of the census machinery):
- positive-cone cusp: t_c = 1510.76411888631, r_c = 3.47459155,
  rho_c = 2.63731321, residual 1.4e-22, G''' = 749.38;
- showcased cusp: t_c = 102.967699396637, residual 1.4e-24,
  G''' = -100.9446 (finite-difference census value -100.9443).

Strengthened unfolding at the certified positive-cone cusp (8 log-spaced
perturbations per ray, both signs, ~2.1 decades):
- rho-ray (generic): exponents 0.340 (+) and 0.324 (-) vs 1/3, with
  POINTWISE prefactor agreement 1-6% against
  (sqrt(3)/2)|6 Z3 drho / G'''|^{1/3} — far cleaner than the first cusp
  (|G'''| = 749 vs 101; better conditioned).
- tangent ray, complex side: exponent 0.495 vs 1/2, prefactor 0.97-1.00;
  real side off = 0 once the three real zeros are resolvable (small-d
  off = 2 readings are line-scan resolution artifacts: three zeros within
  one dt step count as one sign change).
- Incidental structure: near d = +0.016 the excursion collapses to
  5.6e-11 — the tangent ray OSCULATES the true fold curve (the constant
  term is killed only to first order; the quadratic deviation returns the
  ray to the fold locus).  Excluded from the power-law fit as a
  non-generic touching point; itself a normal-form prediction observed.

### Session 7 stratified certification (`src/certify13.py`,
`results/cusp_certification.json`)

- **32/32 stratified cusp candidates certified** by independent 30-digit
  mpmath Newton (residuals <= 6.8e-19).  Census-to-certified displacement:
  median 1.7e-7 in t; two flagged cases both benign on inspection
  (t=897.06: same cusp, census location coarse at ~0.02; t=1975.39:
  displacement 7e-6 RELATIVE to its |rho|=140 coefficients, whence the
  3.4e5 condition number).  |G'''| >= 6.5 across the sample — every
  certified point is a genuine non-degenerate cusp.  Finite-difference
  G''' agrees with mpmath to <= 4.4% (h = 1e-3 stencils).
- **8/8 sample folds** (two per mod-13 family) certified on the 2x2
  system: census fold locations accurate to <= 7e-9 in t.
- Both census-discovered mod-5 zeros verified at |f| ~ 1e-28 and their
  catalog coordinates refined to 30-digit values
  (0.897750154977, 2442.53029246849) and (0.853528651023,
  4392.37642822419).

Defensible summary (tightened in session 7): we observe and numerically
verify, at finite height, the two generic local bifurcation classes
permitted at codimension one and two — folds (sqrt law) and cusps
(cube-root law along generic rays) — with computable normal-form
prefactors confirmed by contour winding at representative points.  NOT
claimed: that every off-line event globally reduces to these classes,
or anything beyond the verified numerical range.
