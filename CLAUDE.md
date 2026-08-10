# DH Collisions — Claude Code briefing

Empirical study of *how zeros leave the critical line* when the Euler product is
absent, using the Davenport–Heilbronn (DH) function as the model organism. Goal:
turn the observed collision mechanism into a precise, testable conjecture — and
connect it to the Weil-positivity picture used in Anthropic's Aug 2026 zeta
lower-bound result (positive/negative-definite subspaces from on/off-line zeros).

## What is established (T ≤ 500, verified; reproduce before extending)

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

## Roadmap (in priority order)

- **Step 1b — the second variable.** Anti-phase + small κ doesn't decide lift-off
  (κ=0.11 stayed on line; κ=0.23 escaped). Hypothesis: the slowly-varying envelope
  of Z_DH at the cancellation point (the same object governing Lehmer pairs)
  decides the branch. Compute it at all 74 anti-phase collisions; test whether
  (κ, envelope) separates outcomes perfectly. Success = a complete local
  mechanism: necessary AND sufficient.
- **Step 2 — statistics at height.** Extend to T ≈ 5000. Use the necessity law to
  make it cheap: enumerate collisions from line scans of the two constituents
  (cheap), then run strip-winding ONLY near anti-phase collisions (~4×+ savings).
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
