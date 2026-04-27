# Missing Inputs

The eval53 large reduction-effect rerun did not start because required
Cantera-readable assets are missing or invalid.

No fallback mechanism, fake trajectory, or copied frozen macro trace was generated.

## Issues

- missing mechanism: benchmarks/assets/eval53_large/mechanisms/diamond_large.yaml
- missing conditions: benchmarks/assets/eval53_large/conditions/diamond.csv
- missing mechanism: benchmarks/assets/eval53_large/mechanisms/sif4_large.yaml
- missing conditions: benchmarks/assets/eval53_large/conditions/sif4.csv
- missing mechanism: benchmarks/assets/eval53_large/mechanisms/ac_large.yaml
- missing conditions: benchmarks/assets/eval53_large/conditions/ac.csv

## Required Layout

- `benchmarks/assets/eval53_large/mechanisms/{benchmark}_large.yaml`
- `benchmarks/assets/eval53_large/conditions/{benchmark}.csv`
- required condition columns: `case_id`, `T`, `P`, `t_end`, `gas_X`

Optional learnCK scores are discovered from:

- `scores/{benchmark}_learnck.csv`
- `learnck_scores/{benchmark}.csv`
- `learnck/{benchmark}_scores.csv`
