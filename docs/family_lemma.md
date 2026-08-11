# The coefficient family G_r = Z1 + r Z2: functional equation lemma

Conventions exactly as implemented (`src/dh_core.py`):

- chi = odd primitive character mod 5 with chi(2) = i; chibar its conjugate.
- tau(chi) = sum_a chi(a) e^{2 pi i a / 5};  eps = tau(chi) / (i sqrt 5),
  |eps| = 1 (odd character); delta = arg(eps)/2 = 0.27678717944852255,
  so eps = e^{2 i delta} and eps_chibar = conj(eps) = e^{-2 i delta}.
- Completed functions Lambda_chi(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) L(s, chi)
  and likewise for chibar.  Functional equations:
      Lambda_chi(s)    = eps      Lambda_chibar(1 - s),
      Lambda_chibar(s) = conj(eps) Lambda_chi(1 - s).
- Coefficient-conjugation identity (Dirichlet series for Re s > 1, then
  analytic continuation, both sides entire):
      conj(Lambda_chi(conj(s))) = Lambda_chibar(s).            (*)

## Lemma 1 (reality of the rotated constituents on the line)

Define on the critical line s = 1/2 + it:

    Ztilde_1(t) = e^{-i delta} Lambda_chi(1/2 + it),
    Ztilde_2(t) = e^{+i delta} Lambda_chibar(1/2 + it).

Both are real for real t.  Proof for Ztilde_1: by (*) and the FE,

    conj(Lambda_chi(1/2 + it)) = Lambda_chibar(1/2 - it)
                               = conj(eps) Lambda_chi(1/2 + it),

so conj(e^{-i delta} Lambda_chi) = e^{i delta} conj(eps) Lambda_chi
= e^{i delta} e^{-2 i delta} Lambda_chi = e^{-i delta} Lambda_chi.  QED.
(Ztilde_2 is the same computation with eps_chibar = e^{-2 i delta}.)

The code's Z_1, Z_2 (`scan2.Z_pair`) are Ztilde_1, Ztilde_2 divided by the
positive factor |(5/pi)^{3/4 + it/2} Gamma(3/4 + it/2)|, i.e.
Z_j(t) = Re[e^{i(theta(t) -/+ delta)} L(1/2+it, chi/chibar)] with
theta(t) = Im log Gamma(3/4 + it/2) + (t/2) log(5/pi).  Positive rescaling
changes no zeros.

## Lemma 2 (Schwarz self-dual FE for every real r)

For real r define

    F_r(s) = e^{-i delta} Lambda_chi(s) + r e^{+i delta} Lambda_chibar(s),

so that F_r(1/2 + it) = Ztilde_1(t) + r Ztilde_2(t) =: G_r(t) (real on the
line by Lemma 1).  Then for ALL complex s:

    F_r(1 - s) = conj(F_r(conj(s))).                           (**)

Proof.  First term: e^{-i delta} Lambda_chi(1-s) = e^{-i delta} eps
Lambda_chibar(s) = e^{i delta} Lambda_chibar(s); and by (*),
conj(e^{-i delta} Lambda_chi(conj(s))) = e^{i delta} Lambda_chibar(s).
Equal.  Second term: r e^{i delta} Lambda_chibar(1-s) = r e^{i delta}
conj(eps) Lambda_chi(s) = r e^{-i delta} Lambda_chi(s); and
conj(r e^{i delta} Lambda_chibar(conj(s))) = r e^{-i delta} Lambda_chi(s)
(r real).  Equal.  QED.

Consequences:
- (**) means rho is a zero iff 1 - conj(rho) is a zero: zeros come in
  mirror pairs (sigma, t) <-> (1 - sigma, t) about the critical LINE at the
  same height.  This is exactly the symmetry the strip-winding bookkeeping
  (total = line + 2 x right-half) requires, for every real r.
- For r != 1 the Dirichlet coefficients e^{-i delta} chi(n) + r e^{i delta}
  chibar(n) are NOT real, so there is no additional (sigma, t) <->
  (sigma, -t) symmetry; all work is at t > 0, as before.  r = 1 is DH.
- r = 0 vertex: F_0 = e^{-i delta} Lambda_chi is a rotated genuine
  Euler-product L-function; its zeros are those of L(s, chi).  Their being
  on the line is GRH for L(s,chi) — a conjecture; the defensible statement
  is that no off-line zeros are detected in the verified numerical range.
- r < 0: same lemma verbatim (r real suffices).  The r < 0 branch is the
  mirror family Z_1 - |r| Z_2, whose birth points are the in-phase
  Wronskian zeros (see bifurcation census).

Numerical verification: `python3 src/family.py` checks Lemma 1 (reality to
~1e-12), Lemma 2 (relative FE residual of (**) to ~1e-11 at random s in the
strip for r in {0.3, 1.0, 2.5}), and the r = 0 vertex (right-half winding of
L_chi alone = 0 in sample windows).
