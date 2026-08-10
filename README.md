# dh-collisions

Empirical study of the mechanism by which zeros leave the critical line for
self-dual combinations of L-functions without an Euler product, via the
Davenport-Heilbronn function. See CLAUDE.md for full state, roadmap, and
validation protocol.

Quick start:
    cd src
    python3 dh_core.py        # validate numerics (must pass first)
    python3 hunt2.py 75.5 88  # must find sigma=0.808517, t=85.699348
    python3 collisions.py     # rebuild constituent zeros + contingency table
    python3 envelope.py       # Step 1b: two-wave model, retreat criterion
    python3 events.py         # event-level summary + survivor catalog

Requires: numpy, scipy, mpmath, matplotlib.
