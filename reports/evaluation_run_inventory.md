# Evaluation Run Inventory

整理日: 2026-04-27

この台帳は、過去に途中停止または中間ファイル欠落が疑われる縮退評価runを、現在のリポジトリに残っている成果物ベースで分類したものです。

## Status Summary

| group | status | what is available | what is missing | use in report |
|---|---|---|---|---|
| `eval53viz` network figures | complete frozen output | 9 method runs x before/after graph, PNG/SVG/DOT/node map | source summary JSON, artifact manifests | main structural comparison |
| `eval53viz` macro diagnostics | complete frozen output | 9 method runs x macro PNG/SVG/CSV, all-condition CSV | independent reduced Cantera trajectories | trace-based density diagnostic only |
| `eval53viz` diagnostic plots | complete derived output | density-ratio summary, condition plot, worst-condition panel, CSV tables | no new physics data | condition-level readability |
| `eval53viz` end-to-end rerun | incomplete / not rerunnable | frozen figures and CSV only | `diamond/sif4/ac` mechanism YAMLs, reduced YAMLs, run artifacts | do not claim independent validation |
| `eval53_large_reduction_effects` rerun | input-gated / not started | rerun entrypoint and missing-input report | `diamond/sif4/ac` mechanism YAMLs and condition CSVs | ready to run when assets are restored |
| `gri30` Cantera rerun | complete auxiliary run | generated reduced YAML, full/reduced trajectories, six case plots, summary CSV | not the `eval53viz` benchmark | workflow demonstration only |

## Eval53viz Reduction Outputs

Each row has all expected static files: before graph DOT/PNG/SVG/node map, after graph DOT/PNG/SVG/node map, and macro PNG/SVG/CSV.

| benchmark | method | network before/after | macro diagnostic | status |
|---|---|---|---|---|
| diamond | baseline | complete | complete | frozen output usable |
| diamond | learnckpp | complete | complete | frozen output usable |
| diamond | pooling | complete | complete | frozen output usable |
| sif4 | baseline | complete | complete | frozen output usable |
| sif4 | learnckpp | complete | complete | frozen output usable |
| sif4 | pooling | complete | complete | frozen output usable |
| ac | baseline | complete | complete | frozen output usable |
| ac | learnckpp | complete | complete | frozen output usable |
| ac | pooling | complete | complete | frozen output usable |

Primary indexes:

- [eval53viz network index](eval53viz_labeled_networks/index.md)
- [eval53viz macro compare index](eval53viz_labeled_networks/macro_compare_index.md)
- [eval53viz condition diagnostics](eval53viz_labeled_networks/diagnostics/index.md)

## Missing Run Inputs

The current checkout does not contain the source files needed to regenerate the `eval53viz` outputs from scratch:

| expected source | present |
|---|---:|
| `reports/diamond_large_eval_summary_eval53viz.json` | no |
| `reports/sif4_large_eval_summary_eval53viz.json` | no |
| `reports/ac_large_eval_summary_eval53viz.json` | no |
| `artifacts/runs/*` for eval53viz | no |
| reduced Cantera mechanism YAML for `diamond/sif4/ac` | no |
| independent reduced Cantera trajectory for `diamond/sif4/ac` | no |

Because of this, the `eval53viz` figures should be treated as a frozen visualization set. They are useful for comparison and explanation, but not enough for an end-to-end rerun or independent reduced-mechanism validation.

## Eval53 Large Rerun Entrypoint

The rerun entrypoint is implemented at:

- [../tools/run_eval53_large_reduction_effects.py](../tools/run_eval53_large_reduction_effects.py)

Default command:

```powershell
py -3.12 tools/run_eval53_large_reduction_effects.py `
  --suite eval53_large `
  --methods learnck,pooling_proxy `
  --levels mild,medium,strong `
  --assets-dir benchmarks/assets/eval53_large `
  --out reports/eval53_large_reduction_effects
```

Current default output:

- [eval53_large_reduction_effects/missing_inputs.md](eval53_large_reduction_effects/missing_inputs.md)

This is the intended current state: the runner does not create fallback mechanisms or copied reduced trajectories while the real large-suite assets are missing.

## Auxiliary Cantera Rerun

The `gri30` rerun is complete and reproducible in the current environment with `py -3.12`, but it is not the same benchmark as `eval53viz`.

| item | path |
|---|---|
| rerun index | [cantera_rerun_gri30_reaction_pruned/index.md](cantera_rerun_gri30_reaction_pruned/index.md) |
| summary table | [cantera_rerun_gri30_reaction_pruned/summary.csv](cantera_rerun_gri30_reaction_pruned/summary.csv) |
| trajectories | [cantera_rerun_gri30_reaction_pruned/trajectories.csv](cantera_rerun_gri30_reaction_pruned/trajectories.csv) |
| generated reduced mechanism | [cantera_rerun_gri30_reaction_pruned/gri30_reaction_pruned_229_of_325.yaml](cantera_rerun_gri30_reaction_pruned/gri30_reaction_pruned_229_of_325.yaml) |
| generator | [../tools/generate_cantera_rerun_comparison.py](../tools/generate_cantera_rerun_comparison.py) |

Additional smoke output for the new eval53 rerun entrypoint:

- [eval53_large_reduction_effects_smoke/output/index.md](eval53_large_reduction_effects_smoke/output/index.md)
- [eval53_large_reduction_effects_smoke/output/summary_all.csv](eval53_large_reduction_effects_smoke/output/summary_all.csv)
- [eval53_large_reduction_effects_surface_smoke/output/index.md](eval53_large_reduction_effects_surface_smoke/output/index.md)
- [eval53_large_reduction_effects_surface_smoke/output/summary_all.csv](eval53_large_reduction_effects_surface_smoke/output/summary_all.csv)

Interpretation:

- This run proves the full/reduced Cantera trajectory comparison workflow works locally.
- The surface smoke uses Cantera's bundled `diamond.yaml` to verify `surface_phase`, coverage, surface ROP, and deposition proxy output.
- It should not be used as evidence that the `diamond/sif4/ac` reduced mechanisms are valid.
- The `eval53viz` report can cite it only as a workflow check.

## Cleanup Decisions

- Replaced the stale mojibake file `reports/eval53viz_labeled_networks/report.md` with a clean index.
- Kept all generated graph, macro, diagnostic, and rerun outputs because they are referenced by `report.md`.
- Did not delete historical generated outputs; the missing source artifacts make deletion riskier than keeping a labeled frozen set.

## Next Required Run

To finish the stopped evaluation properly, regenerate or recover the following and then rerun the comparison:

1. full Cantera mechanism YAML for `diamond`, `sif4`, and `ac`.
2. reduced mechanism YAML for each method: `baseline`, `learnckpp`, `pooling`.
3. condition tables matching the `eval53viz` cases.
4. independent full/reduced Cantera trajectories for T, P, density, major species, surface coverage, and deposition metric.
5. updated report figures from those independent trajectories.
