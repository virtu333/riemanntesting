# dh-collisions

Empirical study of how zeros of self-dual linear combinations of Dirichlet
L-functions leave the critical line, using the Davenport-Heilbronn function
as the entry point and coefficient families G = sum c_j Z_j of rotated
even/odd character waves as the instrument.

Main components (see CLAUDE.md for full state and FINDINGS.md for results):

- **Family lemma** (docs/family_lemma.md): every real combination of the
  rotated waves satisfies one Schwarz self-dual functional equation, so
  zeros come in mirror pairs about the critical line for every coefficient
  choice; Euler-product L-functions sit at the vertices.
- **Wronskian censuses**: double zeros (folds) of two-wave families are
  zeros of W(Z1, Z2) = Z1 Z2' - Z1' Z2 with critical ratio r* = -Z1/Z2;
  triple-zero (cusp) candidates of three-wave families are zeros of the
  3x3 Wronskian det[(Z_i, Z_i', Z_i'')].  These turn coefficient-space
  bifurcations into 1-D scalar scans, reconstruct the off-line zero
  catalog exactly (mod 5, T <= 5000), and found zeros the earlier
  screening pipeline missed.
- **Normal forms**: fold law |sigma-1/2| = |2 Z2/G''|^{1/2} |r-r*|^{1/2}
  (verified <1% near folds); cusp unfoldings verified at representative
  points (cube root on generic rays, square root on the tangent ray).
- **Crossover**: N_off(T, r) birth diagram over three decades in r;
  Psi(lambda) collapse; Selberg-scale variance growth (conditional
  interpretation).

Quick start (mod-5 pipeline):

    cd src
    python3 dh_core.py        # validate numerics (must pass first)
    python3 hunt2.py 75.5 88  # must find sigma=0.808517, t=85.699348
    python3 family.py         # family lemma checks
    python3 scan2.py 0.5 500 T500   # constituent zero scan
    python3 bifurcate.py 0.5 1250 c0  # Wronskian fold census (shard)
    python3 census_analysis.py       # catalog reconstruction + fold stats
    python3 phi_analysis.py          # crossover function

Mod-13 five-wave simplex:

    python3 mod13.py                  # validation (must pass first)
    python3 bifurcate13.py A 0.5 2500 # fold census, family A
    python3 cusp13.py 0.5 2500        # cusp-candidate census (3x3 Wronskian)
    python3 cusp13_verify.py          # unfolding exponents at a cusp

Requires: numpy, scipy, mpmath.  Data in data/ (zero catalogs, scans,
censuses), results in results/ (FINDINGS.md is the log of record).
