# Eval53 Large Reduction Effects

Fresh Cantera full-vs-reduced rerun outputs for baseline, learnCK, and pooling_proxy candidates.

Failed candidates are retained in `summary_all.csv` and their local `status.json`.

## Outputs

- [summary_all.csv](summary_all.csv)
- [effect_matrix.png](effect_matrix.png)
- [effect_matrix.svg](effect_matrix.svg)

## Effect Curves

- [effect_curve_ac_baseline_proxy_dilute_N2.png](effect_curve_ac_baseline_proxy_dilute_N2.png)
- [effect_curve_ac_learnck_style_proxy_dilute_N2.png](effect_curve_ac_learnck_style_proxy_dilute_N2.png)
- [effect_curve_ac_pooling_proxy_dilute_N2.png](effect_curve_ac_pooling_proxy_dilute_N2.png)
- [effect_curve_sif4_baseline_proxy_add_O2.png](effect_curve_sif4_baseline_proxy_add_O2.png)
- [effect_curve_sif4_learnck_style_proxy_add_O2.png](effect_curve_sif4_learnck_style_proxy_add_O2.png)
- [effect_curve_sif4_pooling_proxy_add_O2.png](effect_curve_sif4_pooling_proxy_add_O2.png)

## Candidate Status

| benchmark | method | level | status | directory |
|---|---|---|---|---|
| ac | baseline_proxy | ratio_1000 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_1000](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_1000) |
| ac | baseline_proxy | ratio_900 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_900](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_900) |
| ac | baseline_proxy | ratio_750 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_750](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_750) |
| ac | baseline_proxy | ratio_600 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_600](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_600) |
| ac | baseline_proxy | ratio_500 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_500](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_500) |
| ac | baseline_proxy | ratio_400 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_400](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_400) |
| ac | baseline_proxy | ratio_300 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_300](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_300) |
| ac | baseline | ratio_200 | failed | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_200](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_200) |
| ac | baseline | ratio_100 | failed | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_100](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_100) |
| ac | baseline | ratio_050 | failed | [reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_050](reports/eval53_method_compression_sweep_ac_sif4/ac/baseline/ratio_050) |
| ac | learnck_style_proxy | ratio_1000 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_1000](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_1000) |
| ac | learnck_style_proxy | ratio_900 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_900](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_900) |
| ac | learnck_style_proxy | ratio_750 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_750](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_750) |
| ac | learnck_style_proxy | ratio_600 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_600](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_600) |
| ac | learnck_style_proxy | ratio_500 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_500](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_500) |
| ac | learnck_style_proxy | ratio_400 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_400](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_400) |
| ac | learnck_style_proxy | ratio_300 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_300](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_300) |
| ac | learnck_style_proxy | ratio_200 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_200](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_200) |
| ac | learnck_style_proxy | ratio_100 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_100](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_100) |
| ac | learnck_style_proxy | ratio_050 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_050](reports/eval53_method_compression_sweep_ac_sif4/ac/learnck/ratio_050) |
| ac | pooling_proxy | ratio_1000 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_1000](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_1000) |
| ac | pooling_proxy | ratio_900 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_900](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_900) |
| ac | pooling_proxy | ratio_750 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_750](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_750) |
| ac | pooling_proxy | ratio_600 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_600](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_600) |
| ac | pooling_proxy | ratio_500 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_500](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_500) |
| ac | pooling_proxy | ratio_400 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_400](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_400) |
| ac | pooling_proxy | ratio_300 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_300](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_300) |
| ac | pooling_proxy | ratio_200 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_200](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_200) |
| ac | pooling_proxy | ratio_100 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_100](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_100) |
| ac | pooling_proxy | ratio_050 | success | [reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_050](reports/eval53_method_compression_sweep_ac_sif4/ac/pooling_proxy/ratio_050) |
| sif4 | baseline_proxy | ratio_1000 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_1000](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_1000) |
| sif4 | baseline_proxy | ratio_900 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_900](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_900) |
| sif4 | baseline_proxy | ratio_750 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_750](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_750) |
| sif4 | baseline_proxy | ratio_600 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_600](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_600) |
| sif4 | baseline_proxy | ratio_500 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_500](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_500) |
| sif4 | baseline_proxy | ratio_400 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_400](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_400) |
| sif4 | baseline_proxy | ratio_300 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_300](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_300) |
| sif4 | baseline | ratio_200 | failed | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_200](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_200) |
| sif4 | baseline | ratio_100 | failed | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_100](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_100) |
| sif4 | baseline | ratio_050 | failed | [reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_050](reports/eval53_method_compression_sweep_ac_sif4/sif4/baseline/ratio_050) |
| sif4 | learnck_style_proxy | ratio_1000 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_1000](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_1000) |
| sif4 | learnck_style_proxy | ratio_900 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_900](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_900) |
| sif4 | learnck_style_proxy | ratio_750 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_750](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_750) |
| sif4 | learnck_style_proxy | ratio_600 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_600](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_600) |
| sif4 | learnck_style_proxy | ratio_500 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_500](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_500) |
| sif4 | learnck_style_proxy | ratio_400 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_400](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_400) |
| sif4 | learnck_style_proxy | ratio_300 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_300](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_300) |
| sif4 | learnck_style_proxy | ratio_200 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_200](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_200) |
| sif4 | learnck_style_proxy | ratio_100 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_100](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_100) |
| sif4 | learnck_style_proxy | ratio_050 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_050](reports/eval53_method_compression_sweep_ac_sif4/sif4/learnck/ratio_050) |
| sif4 | pooling_proxy | ratio_1000 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_1000](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_1000) |
| sif4 | pooling_proxy | ratio_900 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_900](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_900) |
| sif4 | pooling_proxy | ratio_750 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_750](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_750) |
| sif4 | pooling_proxy | ratio_600 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_600](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_600) |
| sif4 | pooling_proxy | ratio_500 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_500](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_500) |
| sif4 | pooling_proxy | ratio_400 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_400](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_400) |
| sif4 | pooling_proxy | ratio_300 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_300](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_300) |
| sif4 | pooling_proxy | ratio_200 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_200](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_200) |
| sif4 | pooling_proxy | ratio_100 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_100](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_100) |
| sif4 | pooling_proxy | ratio_050 | success | [reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_050](reports/eval53_method_compression_sweep_ac_sif4/sif4/pooling_proxy/ratio_050) |
