from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, ScalarFormatter


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "eval53_full_time_evolution"

CASES = [
    {
        "benchmark": "diamond",
        "case_id": "surface0",
        "source": "Cantera bundled diamond surface proxy",
        "path": ROOT
        / "reports"
        / "eval53_large_reduction_effects_surface_smoke"
        / "output"
        / "diamond_smoke"
        / "learnck"
        / "mild"
        / "trajectories_full.csv",
    },
    {
        "benchmark": "sif4",
        "case_id": "add_O2",
        "source": "generated eval53 SiF4 proxy, unreduced full trace",
        "path": ROOT
        / "reports"
        / "eval53_large_learnck_sweep_ac_sif4"
        / "sif4"
        / "learnck"
        / "ratio_1000"
        / "trajectories_full.csv",
    },
    {
        "benchmark": "ac",
        "case_id": "dilute_N2",
        "source": "generated eval53 acetylene proxy, unreduced full trace",
        "path": ROOT
        / "reports"
        / "eval53_large_learnck_sweep_ac_sif4"
        / "ac"
        / "learnck"
        / "ratio_1000"
        / "trajectories_full.csv",
    },
]


def _read_rows(path: Path) -> list[dict[str, float | str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    parsed: list[dict[str, float | str]] = []
    for row in rows:
        parsed_row: dict[str, float | str] = {}
        for key, value in row.items():
            if key == "case_id":
                parsed_row[key] = value
                continue
            try:
                parsed_row[key] = float(value)
            except (TypeError, ValueError):
                parsed_row[key] = float("nan")
        parsed.append(parsed_row)
    return parsed


def _array(rows: list[dict[str, float | str]], key: str) -> np.ndarray:
    return np.asarray([float(row.get(key, float("nan"))) for row in rows], dtype=float)


def _series_columns(rows: list[dict[str, float | str]], prefix: str) -> list[str]:
    if not rows:
        return []
    candidates = [key for key in rows[0] if key.startswith(prefix)]
    scored: list[tuple[float, float, str]] = []
    for key in candidates:
        values = _array(rows, key)
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            continue
        scored.append((float(np.nanmax(np.abs(finite))), float(np.nanmax(finite) - np.nanmin(finite)), key))
    scored.sort(reverse=True)
    return [key for _max_value, _span, key in scored[:6]]


def _format_axis(ax: plt.Axes) -> None:
    ax.grid(alpha=0.25)
    ax.xaxis.set_major_locator(MaxNLocator(5))
    x_formatter = ScalarFormatter(useMathText=False)
    x_formatter.set_powerlimits((-2, 2))
    ax.xaxis.set_major_formatter(x_formatter)
    y_formatter = ScalarFormatter(useOffset=False, useMathText=False)
    y_formatter.set_powerlimits((-3, 4))
    ax.yaxis.set_major_formatter(y_formatter)


def _plot_case(meta: dict[str, object], rows: list[dict[str, float | str]]) -> Path:
    benchmark = str(meta["benchmark"])
    case_id = str(meta["case_id"])
    time = _array(rows, "time_s")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    fig.suptitle(f"{benchmark} / {case_id}: unreduced full Cantera time evolution", fontsize=14)

    axes[0, 0].plot(time, _array(rows, "T"), color="#be4b00", linewidth=2)
    axes[0, 0].set_title("Temperature")
    axes[0, 0].set_xlabel("time [s]")
    axes[0, 0].set_ylabel("T [K]")
    _format_axis(axes[0, 0])

    axes[0, 1].plot(time, _array(rows, "density"), color="#0b7285", linewidth=2)
    axes[0, 1].set_title("Density")
    axes[0, 1].set_xlabel("time [s]")
    axes[0, 1].set_ylabel("density [kg/m3]")
    _format_axis(axes[0, 1])

    for key in _series_columns(rows, "X_"):
        axes[1, 0].plot(time, _array(rows, key), linewidth=1.6, label=key[2:])
    axes[1, 0].set_title("Major gas mole fractions")
    axes[1, 0].set_xlabel("time [s]")
    axes[1, 0].set_ylabel("X [-]")
    _format_axis(axes[1, 0])
    axes[1, 0].legend(fontsize=8, loc="best")

    coverage_cols = _series_columns(rows, "coverage_")
    rate_cols = [key for key in ("wdot_abs_sum", "rop_net_abs_sum", "surface_rop_abs_sum", "deposition_proxy") if key in rows[0]]
    if coverage_cols:
        for key in coverage_cols[:5]:
            axes[1, 1].plot(time, _array(rows, key), linewidth=1.6, label=key.replace("coverage_", ""))
        axes[1, 1].set_title("Surface coverages")
        axes[1, 1].set_ylabel("coverage [-]")
    else:
        for key in rate_cols:
            values = _array(rows, key)
            if np.isfinite(values).any():
                axes[1, 1].plot(time, values, linewidth=1.6, label=key)
        axes[1, 1].set_title("Rate activity proxies")
        axes[1, 1].set_ylabel("proxy value")
    axes[1, 1].set_xlabel("time [s]")
    _format_axis(axes[1, 1])
    axes[1, 1].legend(fontsize=8, loc="best")

    png = OUT / f"full_time_evolution_{benchmark}_{case_id}.png"
    svg = OUT / f"full_time_evolution_{benchmark}_{case_id}.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png


def _plot_overview(case_rows: list[tuple[dict[str, object], list[dict[str, float | str]]]]) -> None:
    fig, axes = plt.subplots(len(case_rows), 3, figsize=(14, 9), constrained_layout=True)
    fig.suptitle("Unreduced full Cantera trajectories used before reduction-effect discussion", fontsize=14)
    for row_index, (meta, rows) in enumerate(case_rows):
        benchmark = str(meta["benchmark"])
        case_id = str(meta["case_id"])
        time = _array(rows, "time_s")
        axes[row_index, 0].plot(time, _array(rows, "T"), color="#be4b00", linewidth=2)
        axes[row_index, 0].set_ylabel(f"{benchmark}\nT [K]")
        _format_axis(axes[row_index, 0])

        axes[row_index, 1].plot(time, _array(rows, "density"), color="#0b7285", linewidth=2)
        axes[row_index, 1].set_ylabel("density")
        _format_axis(axes[row_index, 1])

        species_cols = _series_columns(rows, "X_")[:4]
        for key in species_cols:
            axes[row_index, 2].plot(time, _array(rows, key), linewidth=1.4, label=key[2:])
        axes[row_index, 2].set_title(f"{benchmark} / {case_id}" if row_index == 0 else "")
        _format_axis(axes[row_index, 2])
        axes[row_index, 2].legend(fontsize=7, loc="best")
    for ax in axes[-1, :]:
        ax.set_xlabel("time [s]")
    fig.savefig(OUT / "full_time_evolution_overview.png", dpi=180)
    fig.savefig(OUT / "full_time_evolution_overview.svg")
    plt.close(fig)


def _summary_row(meta: dict[str, object], rows: list[dict[str, float | str]]) -> dict[str, object]:
    t = _array(rows, "T")
    density = _array(rows, "density")
    species_bits = []
    for key in _series_columns(rows, "X_")[:4]:
        values = _array(rows, key)
        species_bits.append(f"{key[2:]}={values[-1]:.4g}")
    return {
        "benchmark": meta["benchmark"],
        "case_id": meta["case_id"],
        "source": meta["source"],
        "n_points": len(rows),
        "time_end_s": float(_array(rows, "time_s")[-1]),
        "T_initial_K": float(t[0]),
        "T_final_K": float(t[-1]),
        "T_delta_K": float(t[-1] - t[0]),
        "density_initial": float(density[0]),
        "density_final": float(density[-1]),
        "density_rel_pct": float((density[-1] - density[0]) / density[0] * 100.0) if density[0] else math.nan,
        "major_species_final": "; ".join(species_bits),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    case_rows: list[tuple[dict[str, object], list[dict[str, float | str]]]] = []
    summary: list[dict[str, object]] = []

    for meta in CASES:
        path = Path(meta["path"])
        if not path.exists():
            raise FileNotFoundError(path)
        rows = _read_rows(path)
        if not rows:
            raise RuntimeError(f"empty trajectory: {path}")
        case_rows.append((meta, rows))
        summary.append(_summary_row(meta, rows))
        _plot_case(meta, rows)

    _plot_overview(case_rows)

    summary_path = OUT / "full_time_evolution_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)

    index = [
        "# Eval53 Full Time Evolution Proxy Summary",
        "",
        "This directory collects unreduced full Cantera trajectories before discussing any reduction effect.",
        "",
        "The original eval53 large Cantera mechanisms are not fully recoverable from the current checkout. "
        "The figures are therefore transparent proxy runs: the diamond panel uses the bundled Cantera "
        "diamond surface example, and the ac/sif4 panels use generated eval53 proxy assets from the formal "
        "learnCK compression sweep.",
        "",
        "## Figures",
        "",
        "- [overview](full_time_evolution_overview.png)",
    ]
    for meta in CASES:
        benchmark = meta["benchmark"]
        case_id = meta["case_id"]
        index.append(f"- [{benchmark} / {case_id}](full_time_evolution_{benchmark}_{case_id}.png)")
    index.extend(["", "## Summary", "", f"- [summary CSV]({summary_path.name})", ""])
    (OUT / "index.md").write_text("\n".join(index), encoding="utf-8")


if __name__ == "__main__":
    main()
