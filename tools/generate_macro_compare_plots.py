from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


RUNS = [
    ("diamond", "baseline", "reports/diamond_large_eval_summary_eval53viz.json"),
    ("diamond", "learnckpp", "reports/diamond_large_eval_summary_eval53viz.json"),
    ("diamond", "pooling", "reports/diamond_large_eval_summary_eval53viz.json"),
    ("sif4", "baseline", "reports/sif4_large_eval_summary_eval53viz.json"),
    ("sif4", "learnckpp", "reports/sif4_large_eval_summary_eval53viz.json"),
    ("sif4", "pooling", "reports/sif4_large_eval_summary_eval53viz.json"),
    ("ac", "baseline", "reports/ac_large_eval_summary_eval53viz.json"),
    ("ac", "learnckpp", "reports/ac_large_eval_summary_eval53viz.json"),
    ("ac", "pooling", "reports/ac_large_eval_summary_eval53viz.json"),
]

TRACE_BY_BENCHMARK = {
    "diamond": "artifacts/traces/diamond_benchmarks_diamond_large_trace_eval6.h5",
    "sif4": "artifacts/traces/sif4_benchmark_large_trace.h5",
    "ac": "artifacts/traces/ac_benchmark_large_trace.h5",
}

NODE_RE = re.compile(r"^\s*([A-Za-z0-9_]+)\s+\[(.*)\];\s*$")
ATTR_RE_TEMPLATE = r'{name}="((?:\\.|[^"\\])*)"'
FORMULA_RE = re.compile(r"([A-Z][a-z]?)([0-9.]*)")

ATOMIC_WEIGHTS_G_MOL = {
    "H": 1.00794,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "Ar": 39.948,
    "Si": 28.0855,
    "F": 18.998403,
}

R_UNIVERSAL = 8.31446261815324


def _attr(attrs: str, name: str) -> str:
    match = re.search(ATTR_RE_TEMPLATE.format(name=re.escape(name)), attrs)
    return match.group(1) if match else ""


def _tooltip_field(tooltip: str, name: str) -> str:
    prefix = f"{name}:"
    for part in tooltip.replace("\\n", "\n").splitlines():
        if part.startswith(prefix):
            return part[len(prefix) :].strip()
    return ""


def _label(attrs: str) -> str:
    raw = _attr(attrs, "label")
    raw = raw.replace(r"\"", '"').replace(r"\\", "\\")
    return raw.split(r"\n", 1)[0].strip()


def _members(attrs: str) -> list[str]:
    tooltip = _attr(attrs, "tooltip")
    raw = _tooltip_field(tooltip, "members")
    if not raw:
        raw = _tooltip_field(tooltip, "name")
    if raw:
        return [x.strip() for x in raw.split(",") if x.strip()]
    label = _label(attrs)
    return [label] if label else []


def after_members_from_dot(path: Path) -> set[str]:
    members: set[str] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = NODE_RE.match(line)
        if not match:
            continue
        node_id, attrs = match.groups()
        if node_id in {"graph", "node", "edge"} or node_id.startswith("lane_"):
            continue
        members.update(_members(attrs))
    return members


def is_gas_species(name: str) -> bool:
    if name.endswith("(s)"):
        return False
    if "_bulk" in name or name.endswith("_bulk"):
        return False
    return True


def clean_formula_name(name: str) -> str:
    if name.upper() == "AR":
        return "Ar"
    cleaned = re.sub(r"\([A-Z]\)$", "", name)
    cleaned = re.sub(r"[_-].*$", "", cleaned)
    return cleaned


def molecular_weight_kg_mol(name: str) -> float:
    formula = clean_formula_name(name)
    total = 0.0
    for element, count_raw in FORMULA_RE.findall(formula):
        if element not in ATOMIC_WEIGHTS_G_MOL:
            continue
        count = float(count_raw) if count_raw else 1.0
        total += ATOMIC_WEIGHTS_G_MOL[element] * count
    if total <= 0.0:
        return math.nan
    return total / 1000.0


def mixture_density_kg_m3(
    *,
    x: np.ndarray,
    pressure_pa: np.ndarray,
    temperature_k: np.ndarray,
    species_names: list[str],
    gas_indices: list[int],
) -> np.ndarray:
    mw = np.array([molecular_weight_kg_mol(species_names[i]) for i in gas_indices], dtype=float)
    good = np.isfinite(mw) & (mw > 0)
    if not np.any(good):
        return np.full((x.shape[0],), np.nan, dtype=float)
    gas_x = np.asarray(x[:, gas_indices], dtype=float)[:, good]
    mw = mw[good]
    sums = np.sum(gas_x, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        gas_x_norm = np.divide(gas_x, sums[:, None], out=np.zeros_like(gas_x), where=sums[:, None] > 0)
    mw_mix = gas_x_norm @ mw
    with np.errstate(invalid="ignore", divide="ignore"):
        return pressure_pa * mw_mix / (R_UNIVERSAL * temperature_k)


def _entry_by_mode(summary_path: Path, mode: str) -> dict:
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    for entry in data["entries"]:
        if entry["mode"] == mode:
            return entry
    raise KeyError(f"mode not found in {summary_path}: {mode}")


def _case_metrics(case_group, species_names: list[str], gas_indices_before: list[int], gas_indices_after: list[int]) -> dict:
    x = np.asarray(case_group["X"], dtype=float)
    temp = np.asarray(case_group["temperature"], dtype=float)
    pressure = np.asarray(case_group["pressure"], dtype=float)
    rho_before = mixture_density_kg_m3(
        x=x,
        pressure_pa=pressure,
        temperature_k=temp,
        species_names=species_names,
        gas_indices=gas_indices_before,
    )
    rho_after = mixture_density_kg_m3(
        x=x,
        pressure_pa=pressure,
        temperature_k=temp,
        species_names=species_names,
        gas_indices=gas_indices_after,
    )
    return {
        "T_last_before": float(temp[-1]),
        "T_last_after": float(temp[-1]),
        "T_max_before": float(np.max(temp)),
        "T_max_after": float(np.max(temp)),
        "P_last_before": float(pressure[-1]),
        "P_last_after": float(pressure[-1]),
        "rho_last_before": float(rho_before[-1]),
        "rho_last_after": float(rho_after[-1]),
        "rho_mean_before": float(np.nanmean(rho_before)),
        "rho_mean_after": float(np.nanmean(rho_after)),
    }


def _plot_metric(ax, cases: list[str], before: list[float], after: list[float], title: str, ylabel: str) -> None:
    x = np.arange(len(cases))
    width = 0.38
    ax.bar(x - width / 2, before, width, label="before", color="#4E79A7")
    ax.bar(x + width / 2, after, width, label="after", color="#F28E2B")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(cases, rotation=35, ha="right", fontsize=8)
    ax.grid(axis="y", alpha=0.25)


def write_plot(path: Path, *, title: str, rows: list[dict[str, object]]) -> None:
    cases = [str(row["case"]) for row in rows]
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), constrained_layout=True)
    fig.suptitle(title, fontsize=13)
    plots = [
        ("T_last", "T last", "K"),
        ("T_max", "T max", "K"),
        ("P_last", "P last", "Pa"),
        ("rho_last", "rho last", "kg/m3"),
        ("rho_mean", "rho mean", "kg/m3"),
    ]
    for ax, (prefix, plot_title, ylabel) in zip(axes.flat, plots):
        before = [float(row[f"{prefix}_before"]) for row in rows]
        after = [float(row[f"{prefix}_after"]) for row in rows]
        _plot_metric(ax, cases, before, after, plot_title, ylabel)
    axes.flat[-1].axis("off")
    axes.flat[0].legend(loc="best")
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--out-dir", default="reports/eval53viz_labeled_networks")
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, object]] = []
    plot_rows: list[dict[str, object]] = []
    for benchmark, mode, summary_rel in RUNS:
        entry = _entry_by_mode(root / summary_rel, mode)
        run_id = str(entry["run_id"])
        trace_path = root / TRACE_BY_BENCHMARK[benchmark]
        after_dot = root / "reports" / run_id / "viz" / "network_reaction_edge_detail.dot"
        after_members = after_members_from_dot(after_dot)
        with h5py.File(trace_path, "r") as handle:
            species_names = [str(x) for x in json.loads(str(handle.attrs["species_names"]))]
            name_to_idx = {name: i for i, name in enumerate(species_names)}
            gas_before = [i for i, name in enumerate(species_names) if is_gas_species(name)]
            gas_after = [
                name_to_idx[name]
                for name in species_names
                if name in after_members and is_gas_species(name)
            ]
            if not gas_after:
                gas_after = gas_before
            rows: list[dict[str, object]] = []
            for case in sorted(handle["cases"].keys()):
                metrics = _case_metrics(handle["cases"][case], species_names, gas_before, gas_after)
                row = {
                    "benchmark": benchmark,
                    "mode": mode,
                    "run_id": run_id,
                    "case": str(case),
                    "gas_species_before": len(gas_before),
                    "gas_species_after_basis": len(gas_after),
                    **metrics,
                }
                rows.append(row)
                all_rows.append(row)
        stem = f"{run_id}_macro_compare"
        write_plot(out_dir / stem, title=f"{benchmark} {mode} macro compare", rows=rows)
        plot_rows.append(
            {
                "benchmark": benchmark,
                "mode": mode,
                "run_id": run_id,
                "png": f"{stem}.png",
                "svg": f"{stem}.svg",
                "csv": f"{stem}.csv",
            }
        )
        with (out_dir / f"{stem}.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    with (out_dir / "macro_compare_all.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)

    lines = [
        "# Eval53viz Macro Compare",
        "",
        "Before uses full gas species in the trace. After uses gas species present in the reduced/merged graph, renormalized for ideal-gas density.",
        "Temperature and pressure are copied from the Cantera trace on both sides because reduced macro trajectories were not stored as independent Cantera runs.",
        "",
        "| benchmark | mode | png | svg | csv |",
        "|---|---|---|---|---|",
    ]
    for row in plot_rows:
        lines.append(
            f"| {row['benchmark']} | {row['mode']} | [{row['png']}]({row['png']}) | "
            f"[{row['svg']}]({row['svg']}) | [{row['csv']}]({row['csv']}) |"
        )
    lines.append("")
    lines.append("[all csv](macro_compare_all.csv)")
    lines.append("")
    (out_dir / "macro_compare_index.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {len(plot_rows)} macro comparison plots to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
