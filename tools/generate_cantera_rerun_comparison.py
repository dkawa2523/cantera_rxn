"""Run independent Cantera full-vs-reduced trajectory comparisons.

The script is intentionally standalone so it can be used when the main
artifact store does not yet contain reduced trajectory artifacts. If a reduced
mechanism is not supplied, it can create a reaction-pruned mechanism from the
full mechanism by keeping reactions with the largest initial net rates of
progress across the requested conditions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MAJOR_SPECIES = ("CH4", "O2", "CO2", "H2O", "CO", "H2", "OH")


@dataclass(frozen=True)
class Condition:
    case_id: str
    T0: float
    P0: float
    phi: float
    t_end: float


def _safe_stem(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "case"


def _read_conditions(path: Path, max_cases: int | None) -> list[Condition]:
    rows: list[Condition] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for index, row in enumerate(reader):
            case_id = (row.get("case_id") or row.get("case") or f"case_{index:03d}").strip()
            T0 = float(row.get("T0") or row.get("T") or row.get("temperature"))
            if row.get("P0_atm"):
                P0 = float(row["P0_atm"]) * ct.one_atm
            elif row.get("P_atm"):
                P0 = float(row["P_atm"]) * ct.one_atm
            else:
                P0 = float(row.get("P0") or row.get("P") or ct.one_atm)
            phi = float(row.get("phi") or 1.0)
            t_end = float(row.get("t_end") or row.get("t_final") or 0.02)
            rows.append(Condition(case_id=case_id, T0=T0, P0=P0, phi=phi, t_end=t_end))
            if max_cases is not None and len(rows) >= max_cases:
                break
    if not rows:
        raise ValueError(f"no conditions found in {path}")
    return rows


def _set_gas_state(gas: ct.Solution, condition: Condition, fuel: str, oxidizer: str) -> None:
    gas.TP = condition.T0, condition.P0
    gas.set_equivalence_ratio(condition.phi, fuel=fuel, oxidizer=oxidizer)


def _score_reactions(
    mechanism: Path,
    conditions: Iterable[Condition],
    fuel: str,
    oxidizer: str,
) -> np.ndarray:
    gas = ct.Solution(str(mechanism))
    scores = np.zeros(gas.n_reactions)
    for condition in conditions:
        _set_gas_state(gas, condition, fuel=fuel, oxidizer=oxidizer)
        scores = np.maximum(scores, np.abs(gas.net_rates_of_progress))
    return scores


def _write_pruned_mechanism(
    full_mechanism: Path,
    reduced_mechanism: Path,
    keep_indices: list[int],
) -> None:
    full = ct.Solution(str(full_mechanism))
    reactions = full.reactions()
    kept_reactions = [reactions[index] for index in keep_indices]
    reduced = ct.Solution(
        thermo=full.thermo_model,
        kinetics=full.kinetics_model,
        species=full.species(),
        reactions=kept_reactions,
        name=f"{full.name}_reaction_pruned",
    )
    reduced.write_yaml(str(reduced_mechanism), header=True)


def _prepare_reduced_mechanism(
    full_mechanism: Path,
    reduced_mechanism: Path | None,
    out_dir: Path,
    conditions: list[Condition],
    fuel: str,
    oxidizer: str,
    keep_fraction: float,
) -> tuple[Path, list[dict[str, object]]]:
    if reduced_mechanism is not None:
        return reduced_mechanism, []

    full = ct.Solution(str(full_mechanism))
    keep_count = max(1, min(full.n_reactions, math.ceil(full.n_reactions * keep_fraction)))
    scores = _score_reactions(full_mechanism, conditions, fuel=fuel, oxidizer=oxidizer)
    order = np.argsort(scores)[::-1]
    reactions = full.reactions()
    selected_set = {int(index) for index in order[:keep_count]}
    duplicate_groups: dict[str, set[int]] = {}
    for index, reaction in enumerate(reactions):
        if getattr(reaction, "duplicate", False):
            duplicate_groups.setdefault(reaction.equation, set()).add(index)
    for indices in duplicate_groups.values():
        if selected_set.intersection(indices):
            selected_set.update(indices)
    selected = sorted(selected_set)
    generated = (
        out_dir
        / f"{full_mechanism.stem}_reaction_pruned_{len(selected)}_of_{full.n_reactions}.yaml"
    )
    _write_pruned_mechanism(full_mechanism, generated, selected)

    rows = []
    for index, score in enumerate(scores):
        rows.append(
            {
                "reaction_index": index,
                "kept": index in selected_set,
                "score_abs_initial_rop": float(score),
                "equation": reactions[index].equation,
            }
        )
    return generated, rows


def _time_grid(t_end: float, n_points: int) -> np.ndarray:
    if n_points < 2:
        return np.array([0.0, t_end])
    start = max(t_end * 1.0e-7, 1.0e-12)
    return np.concatenate(([0.0], np.geomspace(start, t_end, n_points - 1)))


def _simulate(
    mechanism: Path,
    condition: Condition,
    fuel: str,
    oxidizer: str,
    species: tuple[str, ...],
    n_points: int,
) -> list[dict[str, float | str]]:
    gas = ct.Solution(str(mechanism))
    _set_gas_state(gas, condition, fuel=fuel, oxidizer=oxidizer)
    reactor = ct.IdealGasConstPressureReactor(gas, energy="on")
    network = ct.ReactorNet([reactor])
    rows: list[dict[str, float | str]] = []
    for target_time in _time_grid(condition.t_end, n_points):
        if target_time > 0.0:
            network.advance(float(target_time))
        thermo = reactor.thermo
        row: dict[str, float | str] = {
            "time_s": float(target_time),
            "T": float(thermo.T),
            "P": float(thermo.P),
            "density": float(thermo.density),
        }
        for name in species:
            row[f"X_{name}"] = float(thermo[name].X[0]) if name in thermo.species_names else 0.0
        rows.append(row)
    return rows


def _trajectory_metric(
    full_rows: list[dict[str, float | str]],
    reduced_rows: list[dict[str, float | str]],
    condition: Condition,
    species: tuple[str, ...],
) -> dict[str, object]:
    full_t = np.array([float(row["T"]) for row in full_rows])
    red_t = np.array([float(row["T"]) for row in reduced_rows])
    full_rho = np.array([float(row["density"]) for row in full_rows])
    red_rho = np.array([float(row["density"]) for row in reduced_rows])
    eps = 1.0e-30

    def ignition_delay(values: np.ndarray) -> float:
        threshold = condition.T0 + 400.0
        hits = np.where(values >= threshold)[0]
        if len(hits) == 0:
            return float("nan")
        return float(full_rows[int(hits[0])]["time_s"])

    metrics: dict[str, object] = {
        "case_id": condition.case_id,
        "T0": condition.T0,
        "P0": condition.P0,
        "phi": condition.phi,
        "t_end": condition.t_end,
        "T_final_full": float(full_t[-1]),
        "T_final_reduced": float(red_t[-1]),
        "T_final_delta": float(red_t[-1] - full_t[-1]),
        "T_final_rel_pct": float((red_t[-1] - full_t[-1]) / max(abs(full_t[-1]), eps) * 100.0),
        "T_max_full": float(np.max(full_t)),
        "T_max_reduced": float(np.max(red_t)),
        "T_max_delta": float(np.max(red_t) - np.max(full_t)),
        "density_final_full": float(full_rho[-1]),
        "density_final_reduced": float(red_rho[-1]),
        "density_final_rel_pct": float(
            (red_rho[-1] - full_rho[-1]) / max(abs(full_rho[-1]), eps) * 100.0
        ),
        "ignition_delay_full_s": ignition_delay(full_t),
        "ignition_delay_reduced_s": ignition_delay(red_t),
    }
    for name in species:
        key = f"X_{name}"
        full_value = float(full_rows[-1][key])
        reduced_value = float(reduced_rows[-1][key])
        metrics[f"{key}_final_full"] = full_value
        metrics[f"{key}_final_reduced"] = reduced_value
        metrics[f"{key}_final_abs_delta"] = reduced_value - full_value
    return metrics


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _save_plot(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png", dpi=180)
    fig.savefig(out_dir / f"{stem}.svg")
    plt.close(fig)


def _plot_case(
    out_dir: Path,
    case_id: str,
    full_rows: list[dict[str, float | str]],
    reduced_rows: list[dict[str, float | str]],
    species: tuple[str, ...],
) -> None:
    time = np.array([float(row["time_s"]) for row in full_rows])
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.0))
    ax = axes[0][0]
    ax.plot(time, [float(row["T"]) for row in full_rows], label="full", linewidth=2.0)
    ax.plot(time, [float(row["T"]) for row in reduced_rows], label="reduced", linestyle="--")
    ax.set_xscale("log")
    ax.set_ylabel("T [K]")
    ax.set_title("Temperature")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    ax = axes[0][1]
    ax.plot(time, [float(row["density"]) for row in full_rows], label="full", linewidth=2.0)
    ax.plot(time, [float(row["density"]) for row in reduced_rows], label="reduced", linestyle="--")
    ax.set_xscale("log")
    ax.set_ylabel("density [kg/m3]")
    ax.set_title("Density")
    ax.grid(alpha=0.25)

    ax = axes[1][0]
    for name in species[:4]:
        ax.plot(time, [float(row[f"X_{name}"]) for row in full_rows], label=f"{name} full")
        ax.plot(
            time,
            [float(row[f"X_{name}"]) for row in reduced_rows],
            linestyle="--",
            label=f"{name} red",
        )
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1.0e-12)
    ax.set_ylabel("mole fraction")
    ax.set_title("Major species")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=7, ncol=2)

    ax = axes[1][1]
    full_t = np.array([float(row["T"]) for row in full_rows])
    red_t = np.array([float(row["T"]) for row in reduced_rows])
    ax.plot(time, red_t - full_t, color="#C44E52")
    ax.axhline(0.0, color="#202124", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_ylabel("T reduced - full [K]")
    ax.set_title("Temperature error")
    ax.grid(alpha=0.25)

    fig.suptitle(f"Independent Cantera rerun: {case_id}", y=1.02)
    _save_plot(fig, out_dir, f"cantera_rerun_{_safe_stem(case_id)}")


def _plot_summary(out_dir: Path, summary_rows: list[dict[str, object]]) -> None:
    labels = [str(row["case_id"]) for row in summary_rows]
    x = np.arange(len(summary_rows))
    width = 0.28
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))

    axes[0].bar(x, [float(row["T_final_rel_pct"]) for row in summary_rows], color="#4E79A7")
    axes[0].axhline(0.0, color="#202124", linewidth=0.8)
    axes[0].set_title("Final T relative error")
    axes[0].set_ylabel("[%]")

    axes[1].bar(x, [float(row["T_max_delta"]) for row in summary_rows], color="#F28E2B")
    axes[1].axhline(0.0, color="#202124", linewidth=0.8)
    axes[1].set_title("Max T delta")
    axes[1].set_ylabel("[K]")

    axes[2].bar(x, [float(row["density_final_rel_pct"]) for row in summary_rows], color="#59A14F")
    axes[2].axhline(0.0, color="#202124", linewidth=0.8)
    axes[2].set_title("Final density relative error")
    axes[2].set_ylabel("[%]")

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.grid(axis="y", alpha=0.25)
    _save_plot(fig, out_dir, "cantera_rerun_summary")


def _write_index(
    out_dir: Path,
    full_mechanism: Path,
    reduced_mechanism: Path,
    generated_reduced: bool,
    conditions: list[Condition],
    species: tuple[str, ...],
    keep_fraction: float,
) -> None:
    lines = [
        "# Independent Cantera Rerun Comparison",
        "",
        "This directory contains full-vs-reduced trajectories generated by fresh Cantera integrations.",
        "",
        f"- full mechanism: `{full_mechanism}`",
        f"- reduced mechanism: `{reduced_mechanism}`",
        f"- reduced mechanism generated in this run: `{generated_reduced}`",
        f"- conditions: {len(conditions)}",
        f"- major species plotted: {', '.join(species)}",
    ]
    if generated_reduced:
        lines.append(
            f"- auto reduction: kept approximately {keep_fraction:.1%} of reactions by initial absolute ROP score"
        )
    lines.extend(
        [
            "",
            "## Figures",
            "",
            "- [summary png](cantera_rerun_summary.png)",
            "- [summary svg](cantera_rerun_summary.svg)",
            "",
        ]
    )
    for condition in conditions:
        stem = f"cantera_rerun_{_safe_stem(condition.case_id)}"
        lines.append(f"- [{condition.case_id} png]({stem}.png) / [svg]({stem}.svg)")
    lines.extend(
        [
            "",
            "## Tables",
            "",
            "- [summary.csv](summary.csv)",
            "- [trajectories.csv](trajectories.csv)",
            "- [metadata.json](metadata.json)",
        ]
    )
    if generated_reduced:
        lines.append("- [reaction_selection.csv](reaction_selection.csv)")
    (out_dir / "index.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    full_mechanism = (ROOT / args.full_mech).resolve()
    reduced_arg = (ROOT / args.reduced_mech).resolve() if args.reduced_mech else None
    conditions_path = (ROOT / args.conditions).resolve()
    out_dir = (ROOT / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    conditions = _read_conditions(conditions_path, args.max_cases)
    species = tuple(name.strip() for name in args.major_species.split(",") if name.strip())
    reduced_mechanism, reaction_rows = _prepare_reduced_mechanism(
        full_mechanism=full_mechanism,
        reduced_mechanism=reduced_arg,
        out_dir=out_dir,
        conditions=conditions,
        fuel=args.fuel,
        oxidizer=args.oxidizer,
        keep_fraction=args.keep_fraction,
    )

    trajectory_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for condition in conditions:
        full_rows = _simulate(
            full_mechanism,
            condition,
            fuel=args.fuel,
            oxidizer=args.oxidizer,
            species=species,
            n_points=args.n_points,
        )
        reduced_rows = _simulate(
            reduced_mechanism,
            condition,
            fuel=args.fuel,
            oxidizer=args.oxidizer,
            species=species,
            n_points=args.n_points,
        )
        for model, rows in (("full", full_rows), ("reduced", reduced_rows)):
            for row in rows:
                trajectory_rows.append(
                    {
                        "case_id": condition.case_id,
                        "model": model,
                        **row,
                    }
                )
        summary_rows.append(_trajectory_metric(full_rows, reduced_rows, condition, species))
        _plot_case(out_dir, condition.case_id, full_rows, reduced_rows, species)

    _write_csv(out_dir / "trajectories.csv", trajectory_rows)
    _write_csv(out_dir / "summary.csv", summary_rows)
    if reaction_rows:
        _write_csv(out_dir / "reaction_selection.csv", reaction_rows)
    _plot_summary(out_dir, summary_rows)

    metadata = {
        "cantera_version": ct.__version__,
        "full_mechanism": str(full_mechanism),
        "reduced_mechanism": str(reduced_mechanism),
        "generated_reduced": reduced_arg is None,
        "conditions": str(conditions_path),
        "fuel": args.fuel,
        "oxidizer": args.oxidizer,
        "keep_fraction": args.keep_fraction if reduced_arg is None else None,
        "n_points": args.n_points,
        "major_species": species,
    }
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _write_index(
        out_dir,
        full_mechanism=full_mechanism,
        reduced_mechanism=reduced_mechanism,
        generated_reduced=reduced_arg is None,
        conditions=conditions,
        species=species,
        keep_fraction=args.keep_fraction,
    )
    print(f"wrote independent Cantera comparison to {out_dir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-mech", default="benchmarks/assets/mechanisms/gri30.yaml")
    parser.add_argument("--reduced-mech", default=None)
    parser.add_argument("--conditions", default="benchmarks/assets/conditions/gri30_small.csv")
    parser.add_argument("--out-dir", default="reports/cantera_rerun_gri30_reaction_pruned")
    parser.add_argument("--fuel", default="CH4:1")
    parser.add_argument("--oxidizer", default="O2:1, N2:3.76")
    parser.add_argument("--keep-fraction", type=float, default=0.7)
    parser.add_argument("--max-cases", type=int, default=6)
    parser.add_argument("--n-points", type=int, default=180)
    parser.add_argument("--major-species", default=",".join(DEFAULT_MAJOR_SPECIES))
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())
