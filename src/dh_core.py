"""Core machinery: vectorized Hurwitz zeta (Euler-Maclaurin), Dirichlet L mod 5,
Davenport-Heilbronn self-dual combination, completed function, FE check.

Validation targets (must reproduce before trusting any downstream result):
  - hurwitz rel err vs mpmath < 1e-12 at (0.5+30i, a=0.2), (1.3+150i, 0.4), (-0.5+300i, 0.6)
  - FE residual |F(s)-F(1-s)|/|F(s)| < 1e-11 at random s in the strip
  - F(1/2+it) real to ~1e-12 relative
  - First off-line zero of f at sigma=0.808517, t=85.699348 (classical value)
"""
import numpy as np
import mpmath as mp

KMAX = 12
BERN = [float(mp.bernoulli(2 * k)) for k in range(1, KMAX + 1)]

def hurwitz(s, a, N=None):
    """Vectorized Hurwitz zeta zeta(s, a) for array of complex s, scalar a in (0,1].
    Euler-Maclaurin with N ~ |Im s| terms. Accurate to ~1e-12 for t up to a few thousand."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    tmax = np.max(np.abs(s.imag)) if s.size else 0.0
    if N is None:
        N = int(max(30, 1.1 * tmax + 20))
    n = np.arange(N)[:, None] + a
    main = np.sum(n ** (-s[None, :]), axis=0)
    Na = N + a
    tail = Na ** (1 - s) / (s - 1) + 0.5 * Na ** (-s)
    corr = np.zeros_like(s)
    poch = np.ones_like(s)
    fact = 1.0
    for k in range(1, KMAX + 1):
        m = 2 * k
        if k == 1:
            poch = s.copy()
        else:
            poch = poch * (s + m - 3) * (s + m - 2)
        fact *= (m - 1) * m
        corr = corr + BERN[k - 1] / fact * poch * Na ** (-s - m + 1)
    return main + tail + corr

# Odd primitive character mod 5 with chi(2) = i
CHI = {1: 1, 2: 1j, 3: -1j, 4: -1}

def L_chi(s, conj=False):
    """L(s, chi) (or chibar) via 5^{-s} sum of Hurwitz zetas."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    out = np.zeros_like(s)
    for a_res, val in CHI.items():
        v = np.conj(val) if conj else val
        out = out + v * hurwitz(s, a_res / 5.0)
    return 5.0 ** (-s) * out

# Root number eps = tau(chi)/(i sqrt(5)); delta = arg(eps)/2 makes coefficients real
TAU = sum(CHI[a] * np.exp(2j * np.pi * a / 5) for a in CHI)
EPS = TAU / (1j * np.sqrt(5))
DELTA = 0.5 * np.angle(EPS)   # = 0.27678717944852255

def f_dh(s):
    """Davenport-Heilbronn: f = e^{-i delta} L(s,chi) + e^{i delta} L(s,chibar).
    Real Dirichlet coefficients 2*Re(e^{-i delta} chi(n)); self-dual completed FE F(s)=F(1-s)."""
    return np.exp(-1j * DELTA) * L_chi(s) + np.exp(1j * DELTA) * L_chi(s, conj=True)

def gamma_np(z):
    """Vectorized complex Gamma via Lanczos (g=7, n=9). Use scipy.special.loggamma
    for log-space work at large |Im z| (this overflows/underflows past t ~ 400)."""
    p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    z = np.asarray(z, dtype=complex)
    refl = z.real < 0.5
    zz = np.where(refl, 1 - z, z)
    x = np.full_like(zz, p[0])
    for i in range(1, 9):
        x = x + p[i] / (zz - 1 + i)
    t = zz + 6.5
    out = np.sqrt(2 * np.pi) * t ** (zz - 0.5) * np.exp(-t) * x
    out = np.where(refl, np.pi / (np.sin(np.pi * z) * out), out)
    return out

def F_completed(s):
    """Lambda(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s). F(s)=F(1-s); real on the line.
    WARNING: over/underflows at large t — use phase-only log-space form (see hunt2.phase_F)."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    return (5 / np.pi) ** ((s + 1) / 2) * gamma_np((s + 1) / 2) * f_dh(s)

if __name__ == "__main__":
    mp.mp.dps = 25
    errs = []
    for (sig, t, a) in [(0.5, 30.0, 0.2), (1.3, 150.0, 0.4), (-0.5, 300.0, 0.6), (2.0, 500.0, 0.8)]:
        ours = hurwitz(np.array([complex(sig, t)]), a)[0]
        ref = complex(mp.zeta(mp.mpc(sig, t), a))
        errs.append(abs(ours - ref) / abs(ref))
    print("hurwitz rel errs:", ["%.2e" % e for e in errs])
    rng = np.random.default_rng(0)
    ss = rng.uniform(-0.5, 1.5, 5) + 1j * rng.uniform(5, 200, 5)
    fe = np.abs(F_completed(ss) - F_completed(1 - ss)) / np.abs(F_completed(ss))
    print("FE rel resid:", ["%.2e" % e for e in fe])
    zz = F_completed(0.5 + 1j * np.array([20.0, 85.0, 333.0]))
    print("imag/real on line:", ["%.2e" % r for r in np.abs(zz.imag / zz.real)])
    print("delta =", DELTA, " eps =", EPS)
