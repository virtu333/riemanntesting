# DH Collisions — Claude Code briefing

Empirical study of *how zeros leave the critical line* when the Euler product is
absent, using the Davenport–Heilbronn (DH) function as the model organism. Goal:
turn the observed collision mechanism into a precise, testable conjecture — and
connect it to the Weil-positivity picture used in Anthropic's Aug 2026 zeta
lower-bound result (positive/negative-definite subspaces from on/off-line zeros).

## What is established (T ≤ 5000; reproduce before extending)

1. **Machinery** (`src/dh_core.py`, `src/hunt2.py`): vectorized Euler–Maclaurin
   Hurwitz zeta (~1e-12 vs mpmath), DH built from the odd character mod 5,
   FE residual ~1e-13. Zero census: **398 zeros in the strip, 378 on the line,
   10 off-line mirror pairs**, zero bookkeeping mismatches.
   First off-line zero: **σ = 0.808517, t = 85.699348** (matches classical value —
   this is the go/no-go validation anchor).
2. **Structural identity**: Z_DH(t) = Z_χ(t) + Z_χ̄(t) exactly — the DH Z-function
   is a sum of two independent real oscillating waves (rotated single-L Z-functions).
3. **Necessity law (the main finding)**: classify mutual-nearest cross-collisions
   between the zeros of Z_χ and Z_χ̄ by the sign of the slope product:
   in-phase (aligned) vs anti-phase (opposed). At T ≤ 500:
   **230 in-phase collisions → 0 off-line zeros; all 10 off-line pairs arise from
   the 74 anti-phase collisions.** Anti-phase is necessary, NOT sufficient.
4. **Excursion correlation**: cancellation depth κ (two-wave residual/(A+B))
   correlates with excursion σ−½ at r ≈ 0.75; raw collision gap at r ≈ 0.83.
   n = 10 is too small to fit a law. An off-line pair ≡ a wrong-sign extremum of
   Z_DH (a "failed Lehmer pair").
5. **Retreat law (Step 1b, `src/envelope.py` + `src/events.py`)**: local
   two-wave model W = A e^{iφ₁} + B e^{iφ₂} fitted from constituent zeros,
   slopes, and extremum heights only (PCHIP; median 4.8% reconstruction error).
   Phase retreat min Θ′ < 0 — impossible for constant amplitudes with common
   carrier; driven by envelope drift (A′B − AB′) — picks **13 of the 56
   anti-phase cancellation events; those 13 contain all 10 off-line events**
   plus 3 survivors (robust to fit window K = 3/4/5; in-phase control 1/230).
   The survivors are Lehmer-like pairs of Z_DH (pair gaps 0.17–0.27 msp).
   The old hypothesis "(κ, envelope) separates perfectly" is **falsified** —
   survivor t=457.54 cancels deeper (0.028) than off-line t=114.18 (0.053).
   Lift-off = the local zero pair going complex (Im t = −(σ−½)); deciding that
   discriminant sign needs model fidelity ≲ κ, so deep events are undecidable
   from ~5%-accurate slow variables. This sensitivity is intrinsic.
6. **Step 2 at T ≤ 5000 (`src/scan2.py`, `src/step2_*.py`)**: 250 off-line
   zeros total (240 in (500,5000], `data/offline_zeros_5000.json`), zero
   bookkeeping mismatches across 309 wound windows. Necessity chain at scale:
   4148 pairs → 1208 anti-phase → 821 events → 286 retreat (35%) → 237 hit
   (83%) + 49 survivors (17%). Blind wall-to-wall audit of 11% of the range:
   **screen completeness 87%** — 3 misses, all near-threshold no-retreat
   events with the largest κ (0.63–0.72) at *constituent-cluster*
   configurations (≥3 zeros of one wave within ~half its mean spacing +
   amplitude collapse) outside the pairwise two-wave model; the largest
   excursions come from this cluster channel. Widen the screen to
   Θ′min/ω < 0.4 in future runs.
   **Excursion law, zero free parameters: σ − ½ = √(2κ)/ω(t)** — median
   actual/predicted 1.017, r = 0.84 (n = 237); free fit exponents (0.533,
   −0.88) vs predicted (0.5, −1). Growth: uniform rate rejected (p ≈ 0.004);
   N_off ~ T·log²T best, T·log T not excluded; anti-event density ∝ log t,
   retreat fraction flat ≈ 0.35, hit rate 0.76 → 0.87.

## Roadmap (in priority order)

- **Step 1b — DONE (see finding 5).** Outcome: retreat criterion is necessary
  and near-sufficient (13 events ⊃ all 10 off-line); perfect slow-variable
  separation is falsified — the branch is a bifurcation-margin discriminant.
  Follow-ups folded into Steps 2 and 4: (a) test retreat necessity
  out-of-sample at height; (b) does the survivor fraction of retreat events
  (3/13 here) stay O(1)?; (c) survivors (DH Lehmer pairs) deserve their own
  census — they are the near-instability points for the Weil-form picture.
- **Step 2 — DONE (see finding 6).** Deliverables exceeded: 250 zeros, the
  √(2κ)/ω law, growth-law discrimination, 49-survivor census, measured
  completeness. New open thread: the **cluster channel** — lift-off at
  configurations of ≥3 same-wave zeros that the pairwise model misses;
  it supplies the largest excursions, so it likely dominates the tail of
  σ−½ at height. Needs its own local model (three-wave / degenerate-zero).
  Deliverables: ~50–100 off-line zeros; fit σ−½ vs κ scaling; decide whether
  N_off(T) grows like T or T·log²T (the collision model predicts a log-power).
  Spot-check in-phase regions to test necessity out-of-sample.
- **Step 3 — degeneration toward Euler products.** Port to even characters mod 13:
  several genuine L-functions share one functional equation, so real combinations
  form a simplex whose vertices satisfy (G)RH. Measure how off-line zero density
  c(f) switches on as you move away from a vertex — rate, exponent, continuity.
  This is the original target question.
- **Step 4 — theory contact.** Express the necessity law in Weil-quadratic-form
  language: anti-phase collisions should correspond to the negative-definite
  directions. If clean, write it up.

## Engineering notes (learned the hard way)

- Never place a winding contour edge at σ = 0.5 + ε — on-line zeros sit on the
  boundary and adaptive refinement explodes. Count the full symmetric strip
  [-1.5, 2.5] with t1 > 0 (keeps Γ poles and trivial zeros outside) and subtract
  the line count.
- All completed-function work at t ≳ 400 must be phase-only in log space
  (scipy loggamma); Λ itself under/overflows.
- Chunk vectorized evals (~1500 pts); the Hurwitz term matrix is (O(t) × N_pts).
- Line scans: dt = 0.01 was collision-free at T ≤ 500; tighten with T (spacing
  shrinks like 1/log t) and always cross-check counts against strip winding.
- Runs are embarrassingly parallel in t-windows; each 12.5-window at t~500 takes
  ~30–60 s single-core.

## Validation protocol for any new range

1. `python3 src/dh_core.py` — all residuals < 1e-11.
2. Rerun a known window (75.5–88.0) — must find the σ=0.808517 zero to 6 decimals.
3. Every window: total = line + 2 × (localized right-half zeros); investigate any
   mismatch (could be a near-line zero with σ−½ < 1e-3 — that would itself be
   interesting; don't paper over it).
4. Contingency table at T ≤ 500 must reproduce 230/0 and 74/10 before extending.

## Data

- `data/offline_zeros_T500.json` — the 10 verified off-line zeros (right half).
- `data/collision_stats_offline_T500.json` — per-zero collision stats
  (gap/mean-spacing, amplitude ratio ρ, κ, excursion).
- `g1.npy`/`g2.npy` (397 constituent zeros each) regenerate in ~2 min via
  `src/collisions.py`.

## Intellectual honesty rules

- This is empirical mathematics: report exact counts, never smooth over
  bookkeeping mismatches, and distinguish "necessary at T ≤ 500" from "necessary".
- Before claiming novelty for any statement, search: Bombieri–Hejhal (zeros of
  linear combinations), Balanzario & Sánchez-Ortiz (DH zero computations),
  Lehmer pairs literature, and the 2023–25 Baluyot–Goldston–Suriajaya–
  Turnage-Butterbaugh papers.
- Numerical "laws" from n=10 are hypotheses. Say so.
