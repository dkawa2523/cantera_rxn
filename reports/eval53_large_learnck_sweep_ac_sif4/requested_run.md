# Requested LearnCK Compression Sweep

Requested cases:

- `ac:dilute_N2`
- `sif4:add_O2`

Requested sweep:

- method: `learnck`
- keep ratios: `1.0`, `0.9`, `0.75`, `0.6`, `0.5`, `0.4`, `0.3`, `0.2`, `0.1`, `0.05`

Command used:

```powershell
py -3.12 tools/generate_eval53_proxy_assets.py
```

```powershell
py -3.12 tools/run_eval53_large_reduction_effects.py `
  --benchmarks ac,sif4 `
  --assets-dir benchmarks/assets/eval53_large `
  --out reports/eval53_large_learnck_sweep_ac_sif4 `
  --methods learnck `
  --case-ids ac:dilute_N2,sif4:add_O2 `
  --sweep-ratios 1.0,0.9,0.75,0.6,0.5,0.4,0.3,0.2,0.1,0.05 `
  --n-points 120 `
  --score-points 60 `
  --max-network-reactions 40 `
  --no-store-rate-arrays
```

Status:

- Completed with generated Cantera-readable proxy assets.
- `ac_large.yaml` is a local GRI30-based acetylene proxy.
- `sif4_large.yaml` is a synthetic constant-cp SiF4/NH3/O2 gas-phase proxy.
- These are not recovered original eval53 large mechanisms.

See:

- [summary_all.csv](summary_all.csv)
- [effect_curve_ac_learnck_style_proxy_dilute_N2.png](effect_curve_ac_learnck_style_proxy_dilute_N2.png)
- [effect_curve_sif4_learnck_style_proxy_add_O2.png](effect_curve_sif4_learnck_style_proxy_add_O2.png)
