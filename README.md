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
    python3 scan2.py 500 2538 s1   # Step 2: fast dual scan (4 shards)
    python3 step2_screen.py   # pairing + retreat screen over (500,5000]
    python3 step2_wind.py 0 4 # winding shard (then: step2_wind.py merge 4)
    python3 step2_blind.py 3000 3100 b3000   # blind completeness audit
    python3 step2_analysis.py # catalog, excursion law, growth stats

Requires: numpy, scipy, mpmath, matplotlib.
