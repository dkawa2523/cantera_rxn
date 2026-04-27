# Eval53viz Output Index

This directory is the frozen visualization output for the `eval53viz` reduction comparison.

Use the repository-level report for interpretation:

- [../../report.md](../../report.md)
- [../../report.html](../../report.html)

Available outputs:

- [network index](index.md): before/after network PNG/SVG/DOT and node maps.
- [macro compare index](macro_compare_index.md): trace-based macro diagnostic PNG/SVG/CSV.
- [macro compare all CSV](macro_compare_all.csv): condition-level density diagnostic table.
- [diagnostic plot index](diagnostics/index.md): additional condition-level diagnostic plots.

Important limitation:

The original `*_eval_summary_eval53viz.json` source summaries, reduced Cantera mechanism YAML files, and independent reduced Cantera trajectories are not present in this checkout. These files are enough to inspect the frozen visualization output, but not enough to rerun the `eval53viz` reduction evaluation end to end.
