# Eval53 Large Proxy Assets

The original eval53 large Cantera mechanisms are not present in this checkout.
These files were generated as transparent Cantera-readable proxies so the
learnCK compression-sweep workflow can run end-to-end.

- `mechanisms/ac_large.yaml`: local GRI30 gas mechanism copy used as an acetylene proxy.
- `mechanisms/sif4_large.yaml`: synthetic constant-cp SiF4/NH3/O2 gas-phase proxy.
- `conditions/ac.csv`: contains `dilute_N2`.
- `conditions/sif4.csv`: contains `add_O2`.

Do not treat these files as recovered original eval53 large mechanisms.
