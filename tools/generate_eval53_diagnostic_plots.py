"""Generate eval53viz density diagnostic plots without optional dependencies.

The figures are diagnostic proxies. They do not replace independent
full-vs-reduced Cantera trajectory validation.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPORT_DIR = ROOT / "reports" / "eval53viz_labeled_networks"
OUT_DIR = REPORT_DIR / "diagnostics"

METHODS = ("baseline", "learnckpp", "pooling")
BENCHMARKS = ("diamond", "sif4", "ac")
COLORS = {"baseline": "#4E79A7", "learnckpp": "#F28E2B", "pooling": "#59A14F"}
MARKERS = {"baseline": "o", "learnckpp": "s", "pooling": "^"}
OK_PCT = 5.0
CAUTION_PCT = 20.0


def _float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def _load_macro_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with (REPORT_DIR / "macro_compare_all.csv").open("r", encoding="utf-8", newline="") as handle:
        for raw in csv.DictReader(handle):
            rho_mean_before = _float(raw, "rho_mean_before")
            rho_mean_after = _float(raw, "rho_mean_after")
            rho_last_before = _float(raw, "rho_last_before")
            rho_last_after = _float(raw, "rho_last_after")
            density_ratio_mean = rho_mean_after / rho_mean_before
            density_ratio_last = rho_last_after / rho_last_before
            gas_before = int(float(raw["gas_species_before"]))
            gas_after = int(float(raw["gas_species_after_basis"]))
            row: dict[str, object] = {
                **raw,
                "gas_species_before": gas_before,
                "gas_species_after_basis": gas_after,
                "T_last_before": _float(raw, "T_last_before"),
                "T_last_after": _float(raw, "T_last_after"),
                "P_last_before": _float(raw, "P_last_before"),
                "P_last_after": _float(raw, "P_last_after"),
                "rho_mean_before": rho_mean_before,
                "rho_mean_after": rho_mean_after,
                "rho_last_before": rho_last_before,
                "rho_last_after": rho_last_after,
                "density_ratio_mean": density_ratio_mean,
                "density_ratio_last": density_ratio_last,
                "density_rel_delta_mean_pct": (density_ratio_mean - 1.0) * 100.0,
                "density_rel_delta_last_pct": (density_ratio_last - 1.0) * 100.0,
                "gas_basis_ratio": gas_after / gas_before,
                "T_last_delta": _float(raw, "T_last_after") - _float(raw, "T_last_before"),
                "P_last_delta": _float(raw, "P_last_after") - _float(raw, "P_last_before"),
            }
            rows.append(row)
    return rows


def _save(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        fig.tight_layout()
    except RuntimeError:
        pass
    fig.savefig(OUT_DIR / f"{name}.png", dpi=180)
    fig.savefig(OUT_DIR / f"{name}.svg")
    plt.close(fig)


def _severity_color(error_pct: float) -> str:
    if error_pct <= OK_PCT:
        return "#2f9e44"
    if error_pct <= CAUTION_PCT:
        return "#f08c00"
    return "#c92a2a"


def _group_rows(rows: Iterable[dict[str, object]]) -> dict[tuple[str, str], list[dict[str, object]]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["benchmark"]), str(row["mode"]))].append(row)
    return grouped


def _density_error(row: dict[str, object]) -> float:
    return abs(float(row["density_rel_delta_mean_pct"]))


def _decision_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    for (benchmark, mode), group in sorted(_group_rows(rows).items()):
        worst = max(group, key=_density_error)
        summary.append(
            {
                "benchmark": benchmark,
                "mode": mode,
                "cases": len(group),
                "gas_species_before": int(group[0]["gas_species_before"]),
                "gas_species_after_basis": int(group[0]["gas_species_after_basis"]),
                "gas_basis_ratio": float(group[0]["gas_basis_ratio"]),
                "mean_density_ratio": float(np.mean([float(r["density_ratio_mean"]) for r in group])),
                "mean_abs_density_error_pct": float(np.mean([_density_error(r) for r in group])),
                "max_abs_density_error_pct": _density_error(worst),
                "worst_case": str(worst["case"]),
                "worst_density_delta_pct": float(worst["density_rel_delta_mean_pct"]),
            }
        )
    return summary


def plot_density_ratio_by_case(rows: list[dict[str, object]]) -> None:
    ordered = sorted(rows, key=lambda r: (str(r["benchmark"]), str(r["case"]), str(r["mode"])))
    labels = [f"{r['benchmark']}\n{r['case']}\n{r['mode']}" for r in ordered]
    colors = [COLORS[str(r["mode"])] for r in ordered]
    fig, ax = plt.subplots(figsize=(15, 5.2))
    ax.bar(range(len(ordered)), [float(r["density_ratio_mean"]) for r in ordered], color=colors)
    ax.axhline(1.0, color="#202124", linewidth=1.0)
    ax.set_ylabel("rho_mean after / before")
    ax.set_title("Per-condition density diagnostic ratio")
    ax.set_xticks(range(len(ordered)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.grid(axis="y", alpha=0.25)
    _save(fig, "eval53viz_density_ratio_by_case")


def plot_density_ratio_summary(rows: list[dict[str, object]]) -> None:
    summary = _decision_summary(rows)
    ordered = [
        item
        for benchmark in BENCHMARKS
        for method in METHODS
        for item in summary
        if item["benchmark"] == benchmark and item["mode"] == method
    ]
    labels = [f"{r['benchmark']}\n{r['mode']}" for r in ordered]
    colors = [COLORS[str(r["mode"])] for r in ordered]
    grouped = _group_rows(rows)
    means = [float(r["mean_density_ratio"]) for r in ordered]
    mins = [min(float(g["density_ratio_mean"]) for g in grouped[(str(r["benchmark"]), str(r["mode"]))]) for r in ordered]
    maxs = [max(float(g["density_ratio_mean"]) for g in grouped[(str(r["benchmark"]), str(r["mode"]))]) for r in ordered]
    yerr = [np.asarray(means) - np.asarray(mins), np.asarray(maxs) - np.asarray(means)]
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    ax.bar(range(len(ordered)), means, yerr=yerr, color=colors, capsize=4)
    ax.axhline(1.0, color="#202124", linewidth=1.0)
    ax.set_ylabel("rho_mean after / before")
    ax.set_title("Mean and range of density diagnostic ratio")
    ax.set_xticks(range(len(ordered)))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.grid(axis="y", alpha=0.25)
    _save(fig, "eval53viz_density_ratio_summary")


def plot_gas_basis_vs_density(rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    for mode in METHODS:
        subset = [row for row in rows if row["mode"] == mode]
        ax.scatter(
            [float(row["gas_basis_ratio"]) for row in subset],
            [float(row["density_ratio_mean"]) for row in subset],
            label=mode,
            s=58,
            alpha=0.82,
            color=COLORS[mode],
            edgecolor="white",
            linewidth=0.7,
        )
    ax.axhline(1.0, color="#202124", linewidth=1.0)
    ax.set_xlabel("gas basis after / before")
    ax.set_ylabel("rho_mean after / before")
    ax.set_title("Gas-basis retention vs density diagnostic ratio")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    for row in sorted(rows, key=_density_error, reverse=True)[:5]:
        ax.annotate(
            f"{row['benchmark']}/{row['mode']}/{row['case']}",
            (float(row["gas_basis_ratio"]), float(row["density_ratio_mean"])),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
        )
    _save(fig, "eval53viz_gas_basis_vs_density_ratio")


def plot_density_diagnostic_dashboard(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    summary = _decision_summary(rows)
    lookup = {(str(row["benchmark"]), str(row["mode"])): row for row in summary}
    matrix = np.asarray(
        [
            [float(lookup[(benchmark, method)]["max_abs_density_error_pct"]) for method in METHODS]
            for benchmark in BENCHMARKS
        ]
    )
    ranked = sorted(rows, key=_density_error, reverse=True)[:12]

    fig = plt.figure(figsize=(14.5, 9.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], width_ratios=[1.0, 1.35])
    ax_heat = fig.add_subplot(grid[0, 0])
    ax_bar = fig.add_subplot(grid[0, 1])
    ax_scatter = fig.add_subplot(grid[1, :])
    fig.suptitle(
        "Eval53viz density diagnostic: which frozen reductions need Cantera rerun first?",
        fontsize=14,
    )

    image = ax_heat.imshow(matrix, cmap="RdYlGn_r", vmin=0.0, vmax=75.0)
    ax_heat.set_title("Worst |rho_after/rho_before - 1| by run")
    ax_heat.set_xticks(range(len(METHODS)))
    ax_heat.set_xticklabels(METHODS, rotation=25, ha="right")
    ax_heat.set_yticks(range(len(BENCHMARKS)))
    ax_heat.set_yticklabels(BENCHMARKS)
    for y, benchmark in enumerate(BENCHMARKS):
        for x, _method in enumerate(METHODS):
            value = matrix[y, x]
            ax_heat.text(
                x,
                y,
                f"{value:.1f}%",
                ha="center",
                va="center",
                color="white" if value > CAUTION_PCT else "#1f2933",
                fontsize=10,
                fontweight="bold",
            )
    cbar = fig.colorbar(image, ax=ax_heat, shrink=0.86)
    cbar.set_label("max absolute density diagnostic error [%]")

    labels = [f"{r['benchmark']}/{r['mode']}\n{r['case']}" for r in ranked]
    errors = [_density_error(row) for row in ranked]
    y_positions = np.arange(len(ranked))
    ax_bar.barh(y_positions, errors, color=[_severity_color(value) for value in errors])
    ax_bar.axvline(OK_PCT, color="#2f9e44", linestyle="--", linewidth=1.0)
    ax_bar.axvline(CAUTION_PCT, color="#f08c00", linestyle="--", linewidth=1.0)
    ax_bar.set_yticks(y_positions)
    ax_bar.set_yticklabels(labels, fontsize=8)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("absolute density diagnostic error [%]")
    ax_bar.set_title("Worst cases, sorted by action priority")
    ax_bar.grid(axis="x", alpha=0.25)
    for y, row in zip(y_positions, ranked):
        value = float(row["density_rel_delta_mean_pct"])
        ax_bar.text(abs(value) + 1.0, y, f"{value:+.1f}%", va="center", fontsize=8)

    max_error = max(errors) if errors else 0.0
    y_max = max(75.0, max_error + 7.0)
    ax_scatter.axhspan(0, OK_PCT, color="#d3f9d8", alpha=0.5)
    ax_scatter.axhspan(OK_PCT, CAUTION_PCT, color="#fff3bf", alpha=0.45)
    ax_scatter.axhspan(CAUTION_PCT, y_max, color="#ffe3e3", alpha=0.4)
    for mode in METHODS:
        subset = [row for row in rows if row["mode"] == mode]
        ax_scatter.scatter(
            [float(row["gas_basis_ratio"]) * 100.0 for row in subset],
            [_density_error(row) for row in subset],
            label=mode,
            marker=MARKERS[mode],
            s=72,
            alpha=0.85,
            color=COLORS[mode],
            edgecolor="white",
            linewidth=0.8,
        )
    ax_scatter.set_xlabel("retained gas basis after / before [%]")
    ax_scatter.set_ylabel("absolute density diagnostic error [%]")
    ax_scatter.set_title("Low retained basis can be safe or unsafe; density error identifies the risky cases")
    ax_scatter.grid(alpha=0.25)
    ax_scatter.set_ylim(0, y_max)
    ax_scatter.legend(ncol=3, frameon=False, loc="upper right")
    ax_scatter.text(99.0, OK_PCT / 2, "OK <=5%", ha="right", va="center", fontsize=8, color="#2b8a3e")
    ax_scatter.text(99.0, (OK_PCT + CAUTION_PCT) / 2, "watch 5-20%", ha="right", va="center", fontsize=8, color="#e67700")
    ax_scatter.text(99.0, CAUTION_PCT + 4, "rerun priority >20%", ha="right", va="center", fontsize=8, color="#c92a2a")

    annotated_keys = {
        ("ac", "learnckpp", "dilute_N2"): "ac/learnckpp/dilute_N2",
        ("sif4", "learnckpp", "add_O2"): "sif4/learnckpp/add_O2",
    }
    for row in rows:
        key = (str(row["benchmark"]), str(row["mode"]), str(row["case"]))
        if key not in annotated_keys:
            continue
        ax_scatter.annotate(
            annotated_keys[key],
            (float(row["gas_basis_ratio"]) * 100.0, _density_error(row)),
            xytext=(6, 8),
            textcoords="offset points",
            fontsize=8,
        )

    _save(fig, "eval53viz_density_diagnostic_dashboard")
    return summary


def plot_single_condition_panels(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    ranked = sorted(rows, key=_density_error, reverse=True)[:6]
    metrics = [
        ("T_last_before", "T_last_after", "T_last [K]"),
        ("P_last_before", "P_last_after", "P_last [Pa]"),
        ("rho_mean_before", "rho_mean_after", "rho_mean [kg/m3]"),
    ]
    fig, axes = plt.subplots(len(ranked), len(metrics), figsize=(10.5, 2.15 * len(ranked)))
    if len(ranked) == 1:
        axes = axes.reshape(1, -1)
    for row_idx, row in enumerate(ranked):
        row_label = f"{row['benchmark']}/{row['mode']}/{row['case']}"
        for col_idx, (before_col, after_col, title) in enumerate(metrics):
            ax = axes[row_idx][col_idx]
            ax.bar(
                ["before", "after"],
                [float(row[before_col]), float(row[after_col])],
                color=["#4E79A7", "#F28E2B"],
            )
            if row_idx == 0:
                ax.set_title(title)
            if col_idx == 0:
                ax.set_ylabel(row_label, fontsize=8)
            ax.grid(axis="y", alpha=0.25)
    fig.suptitle(
        "Worst single-condition macro diagnostic differences "
        "(not independent Cantera reruns)",
        y=1.0,
        fontsize=12,
    )
    _save(fig, "eval53viz_single_condition_macro_panels")
    return ranked


def _write_csv(path: Path, rows: Iterable[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_summary_tables(
    rows: list[dict[str, object]],
    worst: list[dict[str, object]],
    decision_summary: list[dict[str, object]],
) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    columns = [
        "benchmark",
        "mode",
        "case",
        "gas_species_before",
        "gas_species_after_basis",
        "gas_basis_ratio",
        "density_ratio_mean",
        "density_rel_delta_mean_pct",
        "density_ratio_last",
        "density_rel_delta_last_pct",
        "T_last_delta",
        "P_last_delta",
    ]
    ordered = sorted(rows, key=lambda row: (str(row["benchmark"]), str(row["case"]), str(row["mode"])))
    _write_csv(OUT_DIR / "eval53viz_condition_diagnostics.csv", ordered, columns)
    worst_rows = [{**row, "abs_delta": _density_error(row)} for row in worst]
    _write_csv(OUT_DIR / "eval53viz_worst_condition_diagnostics.csv", worst_rows, columns + ["abs_delta"])
    decision_columns = [
        "benchmark",
        "mode",
        "cases",
        "gas_species_before",
        "gas_species_after_basis",
        "gas_basis_ratio",
        "mean_density_ratio",
        "mean_abs_density_error_pct",
        "max_abs_density_error_pct",
        "worst_case",
        "worst_density_delta_pct",
    ]
    _write_csv(OUT_DIR / "eval53viz_density_decision_summary.csv", decision_summary, decision_columns)


def write_index() -> None:
    payload = """# Eval53viz Diagnostic Plots

These figures use `macro_compare_all.csv`. They are diagnostic proxies only:
reduced-side temperature and pressure are copied from the full trace, and
density is recomputed from the reduced gas basis. They are not independent
full-vs-reduced Cantera reruns.

| figure | png | svg |
|---|---|---|
| decision dashboard | [png](eval53viz_density_diagnostic_dashboard.png) | [svg](eval53viz_density_diagnostic_dashboard.svg) |
| density ratio by condition | [png](eval53viz_density_ratio_by_case.png) | [svg](eval53viz_density_ratio_by_case.svg) |
| density ratio summary | [png](eval53viz_density_ratio_summary.png) | [svg](eval53viz_density_ratio_summary.svg) |
| gas basis vs density ratio | [png](eval53viz_gas_basis_vs_density_ratio.png) | [svg](eval53viz_gas_basis_vs_density_ratio.svg) |
| worst single-condition macro panels | [png](eval53viz_single_condition_macro_panels.png) | [svg](eval53viz_single_condition_macro_panels.svg) |

Tables:

- [condition diagnostics](eval53viz_condition_diagnostics.csv)
- [worst condition diagnostics](eval53viz_worst_condition_diagnostics.csv)
- [density decision summary](eval53viz_density_decision_summary.csv)
"""
    (OUT_DIR / "index.md").write_text(payload, encoding="utf-8")


def main() -> None:
    rows = _load_macro_rows()
    plot_density_ratio_by_case(rows)
    plot_density_ratio_summary(rows)
    plot_gas_basis_vs_density(rows)
    decision_summary = plot_density_diagnostic_dashboard(rows)
    worst = plot_single_condition_panels(rows)
    write_summary_tables(rows, worst, decision_summary)
    write_index()
    print(f"wrote diagnostics to {OUT_DIR}")


if __name__ == "__main__":
    main()
