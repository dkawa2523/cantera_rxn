"""Evaluate learnCK and pooling-proxy reduction effects with fresh Cantera runs.

This script intentionally treats the historical eval53viz outputs as frozen
figures. It requires real Cantera-readable assets before it will run the
diamond/sif4/ac large suite. If those assets are missing, it writes an explicit
missing-input report and stops instead of manufacturing fallback results.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

try:  # pragma: no cover - exercised by runtime environment
    import cantera as ct
except Exception:  # pragma: no cover
    ct = None  # type: ignore[assignment]

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BENCHMARKS = ("diamond", "sif4", "ac")
LEVEL_RATIOS = {"mild": 0.75, "medium": 0.50, "strong": 0.25}
METHODS = {"baseline", "learnck", "pooling_proxy"}
REQUIRED_CONDITION_COLUMNS = ("case_id", "T", "P", "t_end", "gas_X")
COMMON_QOI_SPECIES = (
    "CH4",
    "O2",
    "N2",
    "H2",
    "H2O",
    "CO",
    "CO2",
    "C2H2",
    "C2H4",
    "NH3",
    "SiF4",
)


@dataclass(frozen=True)
class Condition:
    case_id: str
    T: float
    P: float
    t_end: float
    gas_X: str
    surface_phase: str | None = None
    surface_coverages: str | None = None
    area: float | None = None
    volume: float | None = None


@dataclass(frozen=True)
class Candidate:
    benchmark: str
    method: str
    method_label: str
    level: str
    keep_ratio: float
    mechanism: Path
    species_before: int
    species_after: int
    reactions_before: int
    reactions_after: int
    selection_rows: list[dict[str, object]]


@dataclass(frozen=True)
class MechanismView:
    gas: Any
    kinetics: Any
    surface: Any | None = None
    adjacent: tuple[Any, ...] = ()
    surface_phase: str | None = None


class InputError(RuntimeError):
    """Raised when required assets are missing or invalid."""


def _rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve())).replace("\\", "/")
    except ValueError:
        return str(path)


def _safe_stem(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "item"


def _split_csv_arg(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def _parse_levels(value: str) -> list[str]:
    levels = _split_csv_arg(value)
    unknown = [level for level in levels if level not in LEVEL_RATIOS]
    if unknown:
        raise InputError(f"unknown levels: {', '.join(unknown)}")
    return levels


def _ratio_label(ratio: float) -> str:
    value = int(round(ratio * 1000.0))
    return f"ratio_{value:03d}"


def _effect_points(levels_arg: str, sweep_ratios_arg: str | None) -> list[tuple[str, float]]:
    if sweep_ratios_arg:
        points: list[tuple[str, float]] = []
        for raw in _split_csv_arg(sweep_ratios_arg):
            try:
                ratio = float(raw)
            except ValueError as exc:
                raise InputError(f"invalid sweep ratio: {raw}") from exc
            if ratio <= 0.0 or ratio > 1.0:
                raise InputError(f"sweep ratios must be in (0, 1]: {raw}")
            points.append((_ratio_label(ratio), ratio))
        if not points:
            raise InputError("sweep-ratios did not contain any values")
        return points
    return [(level, LEVEL_RATIOS[level]) for level in _parse_levels(levels_arg)]


def _parse_methods(value: str) -> list[str]:
    methods = _split_csv_arg(value)
    unknown = [method for method in methods if method not in METHODS]
    if unknown:
        raise InputError(f"unknown methods: {', '.join(unknown)}")
    return methods


def _parse_case_filters(value: str | None) -> dict[str, set[str]] | None:
    if not value:
        return None
    filters: dict[str, set[str]] = {}
    for raw in _split_csv_arg(value):
        if ":" in raw:
            benchmark, case_id = raw.split(":", 1)
            benchmark = benchmark.strip()
            case_id = case_id.strip()
            if not benchmark or not case_id:
                raise InputError(f"invalid case selector: {raw}")
            filters.setdefault(benchmark, set()).add(case_id)
        else:
            filters.setdefault("*", set()).add(raw)
    return filters


def _filter_conditions(
    benchmark: str,
    conditions: list[Condition],
    filters: dict[str, set[str]] | None,
) -> list[Condition]:
    if filters is None:
        return conditions
    wanted = set(filters.get("*", set())) | set(filters.get(benchmark, set()))
    if not wanted:
        return conditions
    filtered = [condition for condition in conditions if condition.case_id in wanted]
    if not filtered:
        available = ", ".join(condition.case_id for condition in conditions)
        raise InputError(
            f"requested case ids not found for {benchmark}: {sorted(wanted)}; "
            f"available: {available}"
        )
    return filtered


def _require_cantera() -> None:
    if ct is None:
        raise InputError(
            "Cantera is not importable in this Python. Use the planned runtime: py -3.12."
        )


def _read_conditions(path: Path, *, max_cases: int | None = None) -> list[Condition]:
    rows: list[Condition] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        missing = [name for name in REQUIRED_CONDITION_COLUMNS if name not in fields]
        if missing:
            raise InputError(f"{_rel(path)} missing required columns: {', '.join(missing)}")
        for index, row in enumerate(reader):
            case_id = (row.get("case_id") or "").strip()
            gas_x = (row.get("gas_X") or "").strip()
            if not case_id:
                raise InputError(f"{_rel(path)} row {index + 2}: empty case_id")
            if not gas_x:
                raise InputError(f"{_rel(path)} row {index + 2}: empty gas_X")
            try:
                T = float(row.get("T") or "")
                P = float(row.get("P") or "")
                t_end = float(row.get("t_end") or "")
            except ValueError as exc:
                raise InputError(f"{_rel(path)} row {index + 2}: invalid T/P/t_end") from exc
            if T <= 0.0 or P <= 0.0 or t_end <= 0.0:
                raise InputError(f"{_rel(path)} row {index + 2}: T/P/t_end must be positive")
            surface_phase = (row.get("surface_phase") or "").strip() or None
            surface_coverages = (row.get("surface_coverages") or "").strip() or None
            area = _optional_float(row.get("area"))
            volume = _optional_float(row.get("volume"))
            rows.append(
                Condition(
                    case_id=case_id,
                    T=T,
                    P=P,
                    t_end=t_end,
                    gas_X=gas_x,
                    surface_phase=surface_phase,
                    surface_coverages=surface_coverages,
                    area=area,
                    volume=volume,
                )
            )
            if max_cases is not None and len(rows) >= max_cases:
                break
    if not rows:
        raise InputError(f"{_rel(path)} has no condition rows")
    return rows


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    return float(value)


def _asset_paths(assets_dir: Path, benchmark: str) -> tuple[Path, Path]:
    mechanism = assets_dir / "mechanisms" / f"{benchmark}_large.yaml"
    conditions = assets_dir / "conditions" / f"{benchmark}.csv"
    return mechanism, conditions


def _score_candidates(assets_dir: Path, benchmark: str) -> list[Path]:
    return [
        assets_dir / "scores" / f"{benchmark}_learnck.csv",
        assets_dir / "learnck_scores" / f"{benchmark}.csv",
        assets_dir / "learnck" / f"{benchmark}_scores.csv",
    ]


def _surface_phase_from_conditions(conditions: list[Condition]) -> str | None:
    phases = {condition.surface_phase for condition in conditions if condition.surface_phase}
    if len(phases) > 1:
        raise InputError(f"multiple surface_phase values are not supported: {sorted(phases)}")
    return next(iter(phases), None)


def _gas_from_surface(surface: Any) -> Any:
    adjacent = getattr(surface, "adjacent", {})
    if "gas" in adjacent:
        return adjacent["gas"]
    for phase in adjacent.values():
        if "gas" in str(getattr(phase, "thermo_model", "")).lower():
            return phase
    raise InputError(f"surface phase {surface.name!r} has no adjacent gas phase")


def _non_gas_adjacent(surface: Any, gas: Any) -> tuple[Any, ...]:
    return tuple(phase for phase in surface.adjacent.values() if phase.name != gas.name)


def _load_mechanism_view(mechanism: Path, surface_phase: str | None = None) -> MechanismView:
    if surface_phase:
        surface = ct.Interface(str(mechanism), surface_phase)  # type: ignore[union-attr]
        gas = _gas_from_surface(surface)
        return MechanismView(
            gas=gas,
            kinetics=surface,
            surface=surface,
            adjacent=_non_gas_adjacent(surface, gas),
            surface_phase=surface_phase,
        )
    gas = ct.Solution(str(mechanism))  # type: ignore[union-attr]
    return MechanismView(gas=gas, kinetics=gas)


def _write_missing_inputs(out_dir: Path, issues: list[str]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Missing Inputs",
        "",
        "The eval53 large reduction-effect rerun did not start because required",
        "Cantera-readable assets are missing or invalid.",
        "",
        "No fallback mechanism, fake trajectory, or copied frozen macro trace was generated.",
        "",
        "## Issues",
        "",
    ]
    lines.extend(f"- {issue}" for issue in issues)
    lines.extend(
        [
            "",
            "## Required Layout",
            "",
            "- `benchmarks/assets/eval53_large/mechanisms/{benchmark}_large.yaml`",
            "- `benchmarks/assets/eval53_large/conditions/{benchmark}.csv`",
            "- required condition columns: `case_id`, `T`, `P`, `t_end`, `gas_X`",
            "",
            "Optional learnCK scores are discovered from:",
            "",
            "- `scores/{benchmark}_learnck.csv`",
            "- `learnck_scores/{benchmark}.csv`",
            "- `learnck/{benchmark}_scores.csv`",
        ]
    )
    (out_dir / "missing_inputs.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out_dir / "index.md").write_text(
        "# Eval53 Large Reduction Effects\n\n"
        "Status: not run because required inputs are missing.\n\n"
        "See [missing_inputs.md](missing_inputs.md).\n",
        encoding="utf-8",
    )


def _validate_assets(
    assets_dir: Path,
    benchmarks: Iterable[str],
    *,
    max_cases: int | None,
) -> dict[str, tuple[Path, Path, list[Condition]]]:
    issues: list[str] = []
    loaded: dict[str, tuple[Path, Path, list[Condition]]] = {}
    for benchmark in benchmarks:
        mechanism, conditions_path = _asset_paths(assets_dir, benchmark)
        if not mechanism.exists():
            issues.append(f"missing mechanism: {_rel(mechanism)}")
        if not conditions_path.exists():
            issues.append(f"missing conditions: {_rel(conditions_path)}")
        conditions: list[Condition] | None = None
        if conditions_path.exists():
            try:
                conditions = _read_conditions(conditions_path, max_cases=max_cases)
                loaded[benchmark] = (mechanism, conditions_path, conditions)
            except InputError as exc:
                issues.append(str(exc))
        if mechanism.exists() and conditions is not None:
            try:
                view = _load_mechanism_view(
                    mechanism,
                    surface_phase=_surface_phase_from_conditions(conditions),
                )
                if view.gas.n_species <= 0:
                    issues.append(f"mechanism has no gas species: {_rel(mechanism)}")
                if view.kinetics.n_reactions <= 0:
                    phase = view.surface_phase or view.gas.name
                    issues.append(f"mechanism phase {phase!r} has no reactions: {_rel(mechanism)}")
            except Exception as exc:
                issues.append(f"cannot load mechanism {_rel(mechanism)}: {exc}")
    if issues:
        raise InputError("\n".join(issues))
    return loaded


def _time_grid(t_end: float, n_points: int) -> np.ndarray:
    if n_points <= 2:
        return np.array([0.0, t_end], dtype=float)
    start = max(t_end * 1.0e-7, 1.0e-12)
    return np.concatenate(([0.0], np.geomspace(start, t_end, n_points - 1)))


def _set_state(gas: Any, condition: Condition) -> None:
    gas.TPX = condition.T, condition.P, condition.gas_X


def _make_reactor(mechanism: Path, condition: Condition) -> tuple[Any, Any, Any | None]:
    if condition.surface_phase:
        surface = ct.Interface(str(mechanism), condition.surface_phase)  # type: ignore[union-attr]
        gas = _gas_from_surface(surface)
    else:
        surface = None
        gas = ct.Solution(str(mechanism))  # type: ignore[union-attr]
    _set_state(gas, condition)
    reactor = ct.IdealGasConstPressureReactor(gas, energy="on")  # type: ignore[union-attr]
    if condition.volume is not None:
        reactor.volume = condition.volume
    if condition.surface_phase:
        if condition.surface_coverages:
            surface.coverages = condition.surface_coverages
        rsurf = ct.ReactorSurface(surface, reactor)  # type: ignore[union-attr]
        if condition.area is not None:
            rsurf.area = condition.area
    return gas, reactor, surface


def _json_float_list(values: Iterable[float]) -> str:
    return json.dumps([float(value) for value in values], separators=(",", ":"))


def _simulate(
    mechanism: Path,
    condition: Condition,
    *,
    major_species: tuple[str, ...],
    n_points: int,
    store_rate_arrays: bool,
) -> list[dict[str, object]]:
    gas, reactor, surface = _make_reactor(mechanism, condition)
    network = ct.ReactorNet([reactor])  # type: ignore[union-attr]
    rows: list[dict[str, object]] = []
    for target_time in _time_grid(condition.t_end, n_points):
        if target_time > network.time:
            network.advance(float(target_time))
        thermo = reactor.thermo
        row: dict[str, object] = {
            "case_id": condition.case_id,
            "time_s": float(network.time),
            "T": float(thermo.T),
            "P": float(thermo.P),
            "density": float(thermo.density),
        }
        for name in major_species:
            row[f"X_{name}"] = (
                float(thermo[name].X[0]) if name in thermo.species_names else 0.0
            )
        try:
            wdot = np.asarray(thermo.net_production_rates, dtype=float)
            row["wdot_abs_sum"] = float(np.sum(np.abs(wdot)))
            if store_rate_arrays:
                row["wdot_json"] = _json_float_list(wdot)
        except Exception:
            row["wdot_abs_sum"] = float("nan")
        try:
            rop = np.asarray(thermo.net_rates_of_progress, dtype=float)
            row["rop_net_abs_sum"] = float(np.sum(np.abs(rop)))
            if store_rate_arrays:
                row["rop_net_json"] = _json_float_list(rop)
        except Exception:
            row["rop_net_abs_sum"] = float("nan")
        if surface is not None:
            for name, value in zip(surface.species_names, surface.coverages):
                row[f"coverage_{name}"] = float(value)
            try:
                srop = np.asarray(surface.net_rates_of_progress, dtype=float)
                row["surface_rop_abs_sum"] = float(np.sum(np.abs(srop)))
                if store_rate_arrays:
                    row["surface_rop_net_json"] = _json_float_list(srop)
                row["deposition_proxy"] = float(np.sum(surface.net_production_rates))
            except Exception:
                row["surface_rop_abs_sum"] = float("nan")
                row["deposition_proxy"] = float("nan")
        rows.append(row)
    return rows


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _read_score_csv(path: Path, reaction_count: int) -> np.ndarray:
    scores = np.zeros(reaction_count, dtype=float)
    seen = False
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or ())
        if not {"reaction_index", "score"}.issubset(fields):
            raise InputError(f"{_rel(path)} must contain reaction_index and score")
        for row in reader:
            idx = int(row["reaction_index"])
            if idx < 0 or idx >= reaction_count:
                raise InputError(f"{_rel(path)} has reaction_index out of range: {idx}")
            scores[idx] = float(row["score"])
            seen = True
    if not seen:
        raise InputError(f"{_rel(path)} has no score rows")
    return scores


def _score_reactions_by_rop(
    mechanism: Path,
    conditions: list[Condition],
    *,
    n_points: int,
    surface_phase: str | None,
) -> np.ndarray:
    view = _load_mechanism_view(mechanism, surface_phase)
    scores = np.zeros(view.kinetics.n_reactions, dtype=float)
    for condition in conditions:
        gas, reactor, _surface = _make_reactor(mechanism, condition)
        network = ct.ReactorNet([reactor])  # type: ignore[union-attr]
        last_time = 0.0
        for target_time in _time_grid(condition.t_end, n_points):
            if target_time > network.time:
                network.advance(float(target_time))
            dt = max(0.0, float(network.time) - last_time)
            last_time = float(network.time)
            try:
                kinetics = _surface if _surface is not None else reactor.thermo
                rop = np.asarray(kinetics.net_rates_of_progress, dtype=float)
            except Exception:
                continue
            if len(rop) == len(scores):
                scores += np.abs(rop) * max(dt, 1.0e-300)
        del gas
    return scores


def _score_species_by_trace(
    full_rows: list[dict[str, object]],
    species_names: list[str],
) -> dict[str, float]:
    scores = {name: 0.0 for name in species_names}
    for row in full_rows:
        for name in species_names:
            key = f"X_{name}"
            value = row.get(key)
            if value is not None:
                scores[name] = max(scores[name], abs(float(value)))
    return scores


def _duplicate_groups(reactions: list[Any]) -> list[set[int]]:
    groups: dict[str, set[int]] = {}
    for index, reaction in enumerate(reactions):
        if getattr(reaction, "duplicate", False):
            groups.setdefault(str(reaction.equation), set()).add(index)
    return [indices for indices in groups.values() if len(indices) > 1]


def _expand_duplicate_groups(selected: set[int], reactions: list[Any]) -> set[int]:
    expanded = set(selected)
    for group in _duplicate_groups(reactions):
        if expanded.intersection(group):
            expanded.update(group)
    return expanded


def _reaction_species(reaction: Any) -> set[str]:
    names = set(reaction.reactants) | set(reaction.products)
    try:
        names.update(reaction.orders)
    except Exception:
        pass
    third_body = getattr(reaction, "third_body", None)
    if third_body is not None:
        try:
            names.update(third_body.efficiencies)
        except Exception:
            pass
    return {str(name) for name in names}


def _write_reduced_solution(
    full_mechanism: Path,
    target: Path,
    *,
    species_names: set[str] | None,
    reaction_indices: set[int],
    name_suffix: str,
    surface_phase: str | None,
) -> tuple[int, int]:
    view = _load_mechanism_view(full_mechanism, surface_phase)
    full = view.kinetics
    reactions = full.reactions()
    reaction_indices = _expand_duplicate_groups(set(reaction_indices), reactions)
    kept_reactions = [reactions[idx] for idx in sorted(reaction_indices)]
    target.parent.mkdir(parents=True, exist_ok=True)
    if view.surface is None:
        if species_names is None:
            species = view.gas.species()
        else:
            needed = set(species_names)
            for idx in reaction_indices:
                needed.update(_reaction_species(reactions[idx]))
            species = [sp for sp in view.gas.species() if sp.name in needed]
        reduced = ct.Solution(  # type: ignore[union-attr]
            thermo=view.gas.thermo_model,
            kinetics=view.gas.kinetics_model,
            species=species,
            reactions=kept_reactions,
            name=f"{view.gas.name}_{name_suffix}",
        )
        reduced.name = f"{view.gas.name}_{name_suffix}"
        reduced.write_yaml(str(target), header=True)
        loaded = ct.Solution(str(target))  # type: ignore[union-attr]
        return int(loaded.n_species), int(loaded.n_reactions)

    gas_names = set(view.gas.species_names) if species_names is None else set(species_names)
    gas_names &= set(view.gas.species_names)
    gas_species = [sp for sp in view.gas.species() if sp.name in gas_names]
    gas_reactions = [
        reaction
        for reaction in view.gas.reactions()
        if _reaction_species(reaction).issubset(gas_names)
    ]
    reduced_gas = ct.Solution(  # type: ignore[union-attr]
        thermo=view.gas.thermo_model,
        kinetics=view.gas.kinetics_model,
        species=gas_species,
        reactions=gas_reactions,
        name=view.gas.name,
    )
    reduced_gas.name = view.gas.name
    adjacent = (reduced_gas, *view.adjacent)
    reduced_surface = ct.Interface(  # type: ignore[union-attr]
        thermo=view.surface.thermo_model,
        kinetics=view.surface.kinetics_model,
        species=view.surface.species(),
        reactions=kept_reactions,
        adjacent=list(adjacent),
        name=view.surface.name,
    )
    reduced_surface.name = view.surface.name
    reduced_surface.write_yaml(str(target), phases=[*adjacent, reduced_surface])
    _inject_adjacent_phases(target, view.surface.name, [phase.name for phase in adjacent])
    loaded_surface = ct.Interface(str(target), view.surface.name)  # type: ignore[union-attr]
    loaded_gas = _gas_from_surface(loaded_surface)
    species_total = int(loaded_gas.n_species + loaded_surface.n_species)
    species_total += sum(int(phase.n_species) for phase in _non_gas_adjacent(loaded_surface, loaded_gas))
    return species_total, int(loaded_surface.n_reactions)


def _learnck_candidate(
    benchmark: str,
    full_mechanism: Path,
    out_dir: Path,
    level: str,
    ratio: float,
    scores: np.ndarray,
    score_source: str,
    surface_phase: str | None,
) -> Candidate:
    view = _load_mechanism_view(full_mechanism, surface_phase)
    full = view.kinetics
    keep_count = max(1, min(full.n_reactions, int(math.ceil(full.n_reactions * ratio))))
    order = np.argsort(np.abs(scores))[::-1]
    selected = {int(idx) for idx in order[:keep_count]}
    mechanism_target = out_dir / "mechanism_reduced.yaml"
    species_after, reactions_after = _write_reduced_solution(
        full_mechanism,
        mechanism_target,
        species_names=None,
        reaction_indices=selected,
        name_suffix=f"learnck_{level}",
        surface_phase=surface_phase,
    )
    selected = _expand_duplicate_groups(selected, full.reactions())
    rows = [
        {
            "kind": "reaction",
            "reaction_index": idx,
            "kept": idx in selected,
            "score": float(scores[idx]),
            "score_source": score_source,
            "equation": full.reaction(idx).equation,
        }
        for idx in range(full.n_reactions)
    ]
    species_before = view.gas.n_species
    if view.surface is not None:
        species_before += view.surface.n_species + sum(phase.n_species for phase in view.adjacent)
    return Candidate(
        benchmark=benchmark,
        method="learnck",
        method_label="learnck" if score_source != "integrated_abs_rop_proxy" else "learnck_style_proxy",
        level=level,
        keep_ratio=ratio,
        mechanism=mechanism_target,
        species_before=species_before,
        species_after=species_after,
        reactions_before=full.n_reactions,
        reactions_after=reactions_after,
        selection_rows=rows,
    )


def _baseline_candidate(
    benchmark: str,
    full_mechanism: Path,
    out_dir: Path,
    level: str,
    ratio: float,
    conditions: list[Condition],
    species_scores: dict[str, float],
    protected_species: set[str],
    surface_phase: str | None,
) -> Candidate:
    view = _load_mechanism_view(full_mechanism, surface_phase)
    full = view.kinetics
    gas_species = list(view.gas.species_names)
    protected = {name for name in protected_species if name in gas_species}
    for condition in conditions:
        gas, _reactor, _surface = _make_reactor(full_mechanism, condition)
        protected.update(name for name, value in zip(gas.species_names, gas.X) if value > 0.0)
    target_species_count = max(len(protected), 1, int(math.ceil(len(gas_species) * ratio)))
    ranked_species = sorted(
        gas_species,
        key=lambda name: (name not in protected, -species_scores.get(name, 0.0), name),
    )
    selected_species = (set(ranked_species[:target_species_count]) | protected) & set(gas_species)
    retained_names = set(selected_species)
    if view.surface is not None:
        retained_names.update(view.surface.species_names)
        for phase in view.adjacent:
            retained_names.update(phase.species_names)
    reactions = full.reactions()
    selected_reactions = {
        idx
        for idx, reaction in enumerate(reactions)
        if _reaction_species(reaction).issubset(retained_names)
    }
    selected_reactions = _expand_duplicate_groups(selected_reactions, reactions)
    if not selected_reactions:
        raise InputError("baseline candidate would remove all reactions")
    mechanism_target = out_dir / "mechanism_reduced.yaml"
    species_after, reactions_after = _write_reduced_solution(
        full_mechanism,
        mechanism_target,
        species_names=selected_species,
        reaction_indices=selected_reactions,
        name_suffix=f"baseline_{level}",
        surface_phase=surface_phase,
    )
    rows: list[dict[str, object]] = []
    for name in gas_species:
        reason = []
        if name in protected:
            reason.append("protected")
        if name in selected_species:
            reason.append("selected")
        rows.append(
            {
                "kind": "species",
                "species": name,
                "kept": name in selected_species,
                "score": float(species_scores.get(name, 0.0)),
                "reason": ",".join(reason) or "removed",
            }
        )
    for idx, reaction in enumerate(reactions):
        rows.append(
            {
                "kind": "reaction",
                "reaction_index": idx,
                "kept": idx in selected_reactions,
                "equation": reaction.equation,
            }
        )
    species_before = view.gas.n_species
    if view.surface is not None:
        species_before += view.surface.n_species + sum(phase.n_species for phase in view.adjacent)
    return Candidate(
        benchmark=benchmark,
        method="baseline",
        method_label="baseline_proxy",
        level=level,
        keep_ratio=ratio,
        mechanism=mechanism_target,
        species_before=species_before,
        species_after=species_after,
        reactions_before=full.n_reactions,
        reactions_after=reactions_after,
        selection_rows=rows,
    )


def _formula_key(species: Any) -> tuple[tuple[str, float], ...]:
    return tuple(sorted((str(k), float(v)) for k, v in species.composition.items()))


def _pooling_candidate(
    benchmark: str,
    full_mechanism: Path,
    out_dir: Path,
    level: str,
    ratio: float,
    conditions: list[Condition],
    species_scores: dict[str, float],
    protected_species: set[str],
    surface_phase: str | None,
) -> Candidate:
    view = _load_mechanism_view(full_mechanism, surface_phase)
    full = view.kinetics
    all_species = view.gas.species()
    species_names = [sp.name for sp in all_species]
    protected = {name for name in protected_species if name in species_names}
    for condition in conditions:
        gas, _reactor, _surface = _make_reactor(full_mechanism, condition)
        protected.update(name for name, value in zip(gas.species_names, gas.X) if value > 0.0)
    target_count = max(len(protected), 1, int(math.ceil(len(species_names) * ratio)))
    ranked = sorted(
        species_names,
        key=lambda name: (name not in protected, -species_scores.get(name, 0.0), name),
    )
    selected_species = set(ranked[:target_count]) | protected
    selected_species &= set(species_names)
    retained_names = set(selected_species)
    if view.surface is not None:
        retained_names.update(view.surface.species_names)
        for phase in view.adjacent:
            retained_names.update(phase.species_names)
    reactions = full.reactions()
    selected_reactions = {
        idx
        for idx, reaction in enumerate(reactions)
        if _reaction_species(reaction).issubset(retained_names)
    }
    selected_reactions = _expand_duplicate_groups(selected_reactions, reactions)
    mechanism_target = out_dir / "mechanism_reduced.yaml"
    species_after, reactions_after = _write_reduced_solution(
        full_mechanism,
        mechanism_target,
        species_names=selected_species,
        reaction_indices=selected_reactions,
        name_suffix=f"pooling_proxy_{level}",
        surface_phase=surface_phase,
    )
    formula_groups: dict[tuple[tuple[str, float], ...], list[str]] = {}
    for sp in all_species:
        formula_groups.setdefault(_formula_key(sp), []).append(sp.name)
    representative_by_group: dict[str, str] = {}
    for members in formula_groups.values():
        representative = max(members, key=lambda name: (species_scores.get(name, 0.0), name))
        for member in members:
            representative_by_group[member] = representative
    rows: list[dict[str, object]] = []
    for name in species_names:
        reason = []
        if name in protected:
            reason.append("protected")
        if name in selected_species:
            reason.append("selected")
        rows.append(
            {
                "kind": "species",
                "species": name,
                "kept": name in selected_species,
                "representative": representative_by_group.get(name, name),
                "score": float(species_scores.get(name, 0.0)),
                "reason": ",".join(reason) or "removed",
            }
        )
    for idx, reaction in enumerate(reactions):
        rows.append(
            {
                "kind": "reaction",
                "reaction_index": idx,
                "kept": idx in selected_reactions,
                "equation": reaction.equation,
            }
        )
    species_before = view.gas.n_species
    if view.surface is not None:
        species_before += view.surface.n_species + sum(phase.n_species for phase in view.adjacent)
    return Candidate(
        benchmark=benchmark,
        method="pooling_proxy",
        method_label="pooling_proxy",
        level=level,
        keep_ratio=ratio,
        mechanism=mechanism_target,
        species_before=species_before,
        species_after=species_after,
        reactions_before=full.n_reactions,
        reactions_after=reactions_after,
        selection_rows=rows,
    )


def _metric_summary(
    benchmark: str,
    candidate: Candidate,
    condition: Condition,
    full_rows: list[dict[str, object]],
    reduced_rows: list[dict[str, object]],
    major_species: tuple[str, ...],
) -> dict[str, object]:
    full_t = np.array([float(row["T"]) for row in full_rows])
    red_t = np.array([float(row["T"]) for row in reduced_rows])
    full_rho = np.array([float(row["density"]) for row in full_rows])
    red_rho = np.array([float(row["density"]) for row in reduced_rows])
    eps = 1.0e-300

    def ignition_delay(rows: list[dict[str, object]], values: np.ndarray) -> float:
        threshold = condition.T + 400.0
        hits = np.where(values >= threshold)[0]
        if len(hits) == 0:
            return float("nan")
        return float(rows[int(hits[0])]["time_s"])

    row: dict[str, object] = {
        "benchmark": benchmark,
        "method": candidate.method,
        "method_label": candidate.method_label,
        "level": candidate.level,
        "keep_ratio_target": candidate.keep_ratio,
        "case_id": condition.case_id,
        "species_before": candidate.species_before,
        "species_after": candidate.species_after,
        "reactions_before": candidate.reactions_before,
        "reactions_after": candidate.reactions_after,
        "species_reduction_pct": float(
            (1.0 - candidate.species_after / max(candidate.species_before, 1)) * 100.0
        ),
        "reaction_reduction_pct": float(
            (1.0 - candidate.reactions_after / max(candidate.reactions_before, 1)) * 100.0
        ),
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
        "ignition_delay_full_s": ignition_delay(full_rows, full_t),
        "ignition_delay_reduced_s": ignition_delay(reduced_rows, red_t),
    }
    for name in major_species:
        key = f"X_{name}"
        if key in full_rows[-1] or key in reduced_rows[-1]:
            row[f"{key}_final_abs_delta"] = float(reduced_rows[-1].get(key, 0.0)) - float(
                full_rows[-1].get(key, 0.0)
            )
    return row


def _save_plot(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png", dpi=180)
    fig.savefig(out_dir / f"{stem}.svg")
    plt.close(fig)


def _plot_timeseries(
    out_dir: Path,
    condition: Condition,
    full_rows: list[dict[str, object]],
    reduced_rows: list[dict[str, object]],
    major_species: tuple[str, ...],
) -> None:
    time = np.array([float(row["time_s"]) for row in full_rows])
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.0))
    axes[0][0].plot(time, [float(row["T"]) for row in full_rows], label="full")
    axes[0][0].plot(time, [float(row["T"]) for row in reduced_rows], "--", label="reduced")
    axes[0][0].set_xscale("log")
    axes[0][0].set_title("Temperature")
    axes[0][0].set_ylabel("K")
    axes[0][0].legend(frameon=False)
    axes[0][0].grid(alpha=0.25)

    axes[0][1].plot(time, [float(row["density"]) for row in full_rows], label="full")
    axes[0][1].plot(time, [float(row["density"]) for row in reduced_rows], "--", label="reduced")
    axes[0][1].set_xscale("log")
    axes[0][1].set_title("Density")
    axes[0][1].set_ylabel("kg/m3")
    axes[0][1].grid(alpha=0.25)

    for name in major_species[:5]:
        key = f"X_{name}"
        axes[1][0].plot(time, [float(row.get(key, 0.0)) for row in full_rows], label=f"{name} full")
        axes[1][0].plot(
            time,
            [float(row.get(key, 0.0)) for row in reduced_rows],
            "--",
            label=f"{name} red",
        )
    axes[1][0].set_xscale("log")
    axes[1][0].set_yscale("symlog", linthresh=1.0e-12)
    axes[1][0].set_title("Major species")
    axes[1][0].set_ylabel("mole fraction")
    axes[1][0].legend(frameon=False, fontsize=7, ncol=2)
    axes[1][0].grid(alpha=0.25)

    full_t = np.array([float(row["T"]) for row in full_rows])
    red_t = np.array([float(row["T"]) for row in reduced_rows])
    axes[1][1].plot(time, red_t - full_t, color="#C44E52")
    axes[1][1].axhline(0.0, color="#202124", linewidth=0.8)
    axes[1][1].set_xscale("log")
    axes[1][1].set_title("Temperature error")
    axes[1][1].set_ylabel("K")
    axes[1][1].grid(alpha=0.25)
    _save_plot(fig, out_dir, f"timeseries_{_safe_stem(condition.case_id)}")


def _plot_summary(out_dir: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    labels = [str(row["case_id"]) for row in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))
    axes[0].bar(x, [float(row["T_final_rel_pct"]) for row in rows], color="#4E79A7")
    axes[0].set_title("Final T rel error [%]")
    axes[1].bar(x, [float(row["T_max_delta"]) for row in rows], color="#F28E2B")
    axes[1].set_title("Max T delta [K]")
    axes[2].bar(x, [float(row["density_final_rel_pct"]) for row in rows], color="#59A14F")
    axes[2].set_title("Final density rel error [%]")
    for ax in axes:
        ax.axhline(0.0, color="#202124", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.grid(axis="y", alpha=0.25)
    _save_plot(fig, out_dir, "summary")


def _network_edges(kinetics: Any, max_reactions: int) -> tuple[list[str], list[tuple[str, str, str]]]:
    reactions = kinetics.reactions()
    ranked = sorted(
        range(kinetics.n_reactions),
        key=lambda idx: (len(_reaction_species(reactions[idx])), -idx),
        reverse=True,
    )[:max_reactions]
    species: set[str] = set()
    edges: list[tuple[str, str, str]] = []
    for idx in sorted(ranked):
        reaction = reactions[idx]
        rxn_node = f"R{idx}"
        for name in reaction.reactants:
            species.add(name)
            edges.append((name, rxn_node, "reactant"))
        for name in reaction.products:
            species.add(name)
            edges.append((rxn_node, name, "product"))
    return sorted(species) + [f"R{idx}" for idx in sorted(ranked)], edges


def _write_network(
    mechanism: Path,
    out_dir: Path,
    stem: str,
    title: str,
    *,
    selection_rows: list[dict[str, object]] | None = None,
    max_reactions: int = 40,
    surface_phase: str | None = None,
) -> None:
    view = _load_mechanism_view(mechanism, surface_phase)
    nodes, edges = _network_edges(view.kinetics, max_reactions=max_reactions)
    node_set = set(nodes)
    dot_lines = ["digraph G {", "  rankdir=LR;", f'  label="{title}";']
    for node in nodes:
        shape = "box" if node.startswith("R") else "ellipse"
        dot_lines.append(f'  "{node}" [shape={shape}];')
    for source, target, role in edges:
        if source in node_set and target in node_set:
            dot_lines.append(f'  "{source}" -> "{target}" [label="{role}"];')
    dot_lines.append("}")
    (out_dir / f"{stem}.dot").write_text("\n".join(dot_lines) + "\n", encoding="utf-8")

    species_nodes = [node for node in nodes if not node.startswith("R")]
    reaction_nodes = [node for node in nodes if node.startswith("R")]
    y_species = np.linspace(0.0, 1.0, max(len(species_nodes), 1))
    y_reactions = np.linspace(0.0, 1.0, max(len(reaction_nodes), 1))
    coords = {node: (0.0, float(y_species[i])) for i, node in enumerate(species_nodes)}
    coords.update({node: (1.0, float(y_reactions[i])) for i, node in enumerate(reaction_nodes)})
    fig, ax = plt.subplots(figsize=(9.0, 7.0))
    for source, target, _role in edges:
        if source in coords and target in coords:
            x0, y0 = coords[source]
            x1, y1 = coords[target]
            ax.plot([x0, x1], [y0, y1], color="#8A9199", linewidth=0.5, alpha=0.6)
    for node in species_nodes:
        x, y = coords[node]
        ax.scatter([x], [y], s=42, color="#4E79A7", zorder=3)
        ax.text(x - 0.02, y, node, ha="right", va="center", fontsize=6)
    for node in reaction_nodes:
        x, y = coords[node]
        ax.scatter([x], [y], s=36, color="#F28E2B", marker="s", zorder=3)
        ax.text(x + 0.02, y, node, ha="left", va="center", fontsize=6)
    ax.set_title(title)
    ax.set_axis_off()
    _save_plot(fig, out_dir, stem)

    map_lines = [f"# {title}", "", f"- mechanism: `{_rel(mechanism)}`", ""]
    map_lines.append("## Drawn Nodes")
    map_lines.extend(f"- `{node}`" for node in nodes)
    if selection_rows:
        map_lines.extend(["", "## Selection Map", ""])
        for row in selection_rows:
            if row.get("kind") == "species":
                map_lines.append(
                    f"- {row.get('species')}: kept={row.get('kept')} "
                    f"representative={row.get('representative')} reason={row.get('reason')}"
                )
    (out_dir / f"node_map_{stem.split('_')[-1]}.md").write_text(
        "\n".join(map_lines) + "\n",
        encoding="utf-8",
    )


def _write_status(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _inject_adjacent_phases(path: Path, surface_name: str, adjacent_names: Iterable[str]) -> None:
    adjacent = json.dumps(list(adjacent_names))
    lines = path.read_text(encoding="utf-8").splitlines()
    out: list[str] = []
    in_surface = False
    injected = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("- name: "):
            in_surface = stripped == f"- name: {surface_name}"
        out.append(line)
        if in_surface and stripped.startswith("adjacent-phases:"):
            injected = True
        elif in_surface and stripped.startswith("thermo:") and not injected:
            out.append(f"    adjacent-phases: {adjacent}")
            injected = True
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def _major_species(gas: Any, conditions: list[Condition], requested: str) -> tuple[str, ...]:
    if requested != "auto":
        return tuple(name for name in _split_csv_arg(requested) if name in gas.species_names)
    names: list[str] = []
    for name in COMMON_QOI_SPECIES:
        if name in gas.species_names and name not in names:
            names.append(name)
    for condition in conditions:
        gas.TPX = condition.T, condition.P, condition.gas_X
        for name, value in sorted(zip(gas.species_names, gas.X), key=lambda item: item[1], reverse=True):
            if value > 0.0 and name not in names:
                names.append(name)
            if len(names) >= 10:
                break
    return tuple(names[:10])


def _plot_effect_matrix(out_dir: Path, rows: list[dict[str, object]]) -> None:
    success = [row for row in rows if row.get("status") == "success"]
    if not success:
        return
    groups: dict[tuple[str, str, str], list[float]] = {}
    for row in success:
        key = (str(row["benchmark"]), str(row["method_label"]), str(row["level"]))
        groups.setdefault(key, []).append(abs(float(row.get("T_final_rel_pct", 0.0))))
    labels = ["/".join(key) for key in sorted(groups)]
    values = [float(np.mean(groups[key])) for key in sorted(groups)]
    fig, ax = plt.subplots(figsize=(max(7.0, len(labels) * 0.7), 4.0))
    ax.bar(range(len(labels)), values, color="#4E79A7")
    ax.set_ylabel("mean |final T rel error| [%]")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.25)
    _save_plot(fig, out_dir, "effect_matrix")


def _plot_effect_curves(out_dir: Path, rows: list[dict[str, object]]) -> None:
    success = [row for row in rows if row.get("status") == "success"]
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = {}
    for row in success:
        case_id = str(row.get("case_id", ""))
        if not case_id:
            continue
        key = (
            str(row.get("benchmark")),
            str(row.get("method_label", row.get("method"))),
            case_id,
        )
        groups.setdefault(key, []).append(row)
    for (benchmark, method, case_id), group_rows in groups.items():
        if len(group_rows) < 2:
            continue
        group_rows = sorted(
            group_rows,
            key=lambda row: (
                float(row.get("keep_ratio_target", 0.0)),
                float(row.get("reactions_after", 0.0)),
            ),
        )
        x = np.array([float(row.get("reactions_after", 0.0)) for row in group_rows])
        labels = [str(row.get("level")) for row in group_rows]
        fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2))
        axes[0][0].plot(x, [abs(float(row.get("T_final_rel_pct", 0.0))) for row in group_rows], marker="o")
        axes[0][0].set_ylabel("|final T rel error| [%]")
        axes[0][1].plot(x, [abs(float(row.get("density_final_rel_pct", 0.0))) for row in group_rows], marker="o")
        axes[0][1].set_ylabel("|final density rel error| [%]")
        axes[1][0].plot(x, [float(row.get("reaction_reduction_pct", 0.0)) for row in group_rows], marker="o")
        axes[1][0].set_ylabel("reaction reduction [%]")
        species_keys = [
            key
            for key in group_rows[0]
            if key.startswith("X_") and key.endswith("_final_abs_delta")
        ][:5]
        for key in species_keys:
            axes[1][1].plot(
                x,
                [abs(float(row.get(key, 0.0))) for row in group_rows],
                marker="o",
                label=key.replace("X_", "").replace("_final_abs_delta", ""),
            )
        axes[1][1].set_ylabel("major X final abs error")
        if species_keys:
            axes[1][1].legend(frameon=False, fontsize=7)
        for ax in axes.ravel():
            ax.set_xlabel("kept reactions")
            ax.grid(alpha=0.25)
            if ax.lines:
                ydata = ax.lines[0].get_ydata()
                for idx, (xi, label) in enumerate(zip(x, labels)):
                    ax.annotate(label, (xi, ydata[idx]), fontsize=6)
        fig.suptitle(f"{benchmark} / {method} / {case_id} compression sweep")
        _save_plot(fig, out_dir, f"effect_curve_{_safe_stem(benchmark)}_{_safe_stem(method)}_{_safe_stem(case_id)}")


def _write_index(out_dir: Path, rows: list[dict[str, object]], missing: bool = False) -> None:
    if missing:
        return
    lines = [
        "# Eval53 Large Reduction Effects",
        "",
        "Fresh Cantera full-vs-reduced rerun outputs for baseline, learnCK, and pooling_proxy candidates.",
        "",
        "Failed candidates are retained in `summary_all.csv` and their local `status.json`.",
        "",
        "## Outputs",
        "",
        "- [summary_all.csv](summary_all.csv)",
        "- [effect_matrix.png](effect_matrix.png)",
        "- [effect_matrix.svg](effect_matrix.svg)",
        "",
    ]
    curves = sorted(path.name for path in out_dir.glob("effect_curve_*.png"))
    if curves:
        lines.extend(["## Effect Curves", ""])
        lines.extend(f"- [{name}]({name})" for name in curves)
        lines.append("")
    lines.extend(
        [
            "## Candidate Status",
            "",
            "| benchmark | method | level | status | directory |",
            "|---|---|---|---|---|",
        ]
    )
    seen: set[tuple[str, str, str, str, str]] = set()
    for row in rows:
        key = (
            str(row.get("benchmark")),
            str(row.get("method_label", row.get("method"))),
            str(row.get("level")),
            str(row.get("status")),
            str(row.get("directory", "")),
        )
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | [{key[4]}]({key[4]}) |")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "index.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run_candidate(
    benchmark: str,
    full_mechanism: Path,
    conditions: list[Condition],
    candidate: Candidate,
    candidate_dir: Path,
    *,
    full_cache: dict[str, list[dict[str, object]]],
    major_species: tuple[str, ...],
    n_points: int,
    store_rate_arrays: bool,
    max_network_reactions: int,
    surface_phase: str | None,
) -> list[dict[str, object]]:
    _write_csv(candidate_dir / "selection.csv", candidate.selection_rows)
    _write_network(
        full_mechanism,
        candidate_dir,
        "network_before",
        f"{benchmark} full network",
        max_reactions=max_network_reactions,
        surface_phase=surface_phase,
    )
    _write_network(
        candidate.mechanism,
        candidate_dir,
        "network_after",
        f"{benchmark} {candidate.method_label} {candidate.level} network",
        selection_rows=candidate.selection_rows,
        max_reactions=max_network_reactions,
        surface_phase=surface_phase,
    )
    reduced_rows_all: list[dict[str, object]] = []
    full_rows_all: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for condition in conditions:
        full_rows = full_cache[condition.case_id]
        reduced_rows = _simulate(
            candidate.mechanism,
            condition,
            major_species=major_species,
            n_points=n_points,
            store_rate_arrays=store_rate_arrays,
        )
        full_rows_all.extend(full_rows)
        reduced_rows_all.extend(reduced_rows)
        _plot_timeseries(candidate_dir, condition, full_rows, reduced_rows, major_species)
        summary_rows.append(
            _metric_summary(benchmark, candidate, condition, full_rows, reduced_rows, major_species)
        )
    _write_csv(candidate_dir / "trajectories_full.csv", full_rows_all)
    _write_csv(candidate_dir / "trajectories_reduced.csv", reduced_rows_all)
    _write_csv(candidate_dir / "summary.csv", summary_rows)
    _plot_summary(candidate_dir, summary_rows)
    status = {
        "status": "success",
        "benchmark": benchmark,
        "method": candidate.method,
        "method_label": candidate.method_label,
        "level": candidate.level,
        "mechanism_reduced": _rel(candidate.mechanism),
        "species_before": candidate.species_before,
        "species_after": candidate.species_after,
        "reactions_before": candidate.reactions_before,
        "reactions_after": candidate.reactions_after,
    }
    _write_status(candidate_dir / "status.json", status)
    for row in summary_rows:
        row["status"] = "success"
        row["directory"] = _rel(candidate_dir)
    return summary_rows


def run(args: argparse.Namespace) -> int:
    _require_cantera()
    assets_dir = (ROOT / args.assets_dir).resolve()
    out_dir = (ROOT / args.out).resolve()
    benchmarks = tuple(_split_csv_arg(args.benchmarks))
    methods = _parse_methods(args.methods)
    effect_points = _effect_points(args.levels, args.sweep_ratios)
    case_filters = _parse_case_filters(args.case_ids)
    protected_species = set(_split_csv_arg(args.protected_species))
    try:
        loaded = _validate_assets(assets_dir, benchmarks, max_cases=args.max_cases)
    except InputError as exc:
        issues = [line for line in str(exc).splitlines() if line.strip()]
        _write_missing_inputs(out_dir, issues)
        print(f"missing inputs; wrote {_rel(out_dir / 'missing_inputs.md')}")
        return 2
    try:
        loaded = {
            benchmark: (
                mechanism,
                conditions_path,
                _filter_conditions(benchmark, conditions, case_filters),
            )
            for benchmark, (mechanism, conditions_path, conditions) in loaded.items()
        }
    except InputError as exc:
        _write_missing_inputs(out_dir, [str(exc)])
        print(f"missing inputs; wrote {_rel(out_dir / 'missing_inputs.md')}")
        return 2

    all_rows: list[dict[str, object]] = []
    for benchmark in benchmarks:
        full_mechanism, _conditions_path, conditions = loaded[benchmark]
        surface_phase = _surface_phase_from_conditions(conditions)
        view = _load_mechanism_view(full_mechanism, surface_phase)
        full = view.kinetics
        major_species = _major_species(view.gas, conditions, args.major_species)
        full_cache: dict[str, list[dict[str, object]]] = {}
        full_rows_all: list[dict[str, object]] = []
        for condition in conditions:
            rows = _simulate(
                full_mechanism,
                condition,
                major_species=major_species,
                n_points=args.n_points,
                store_rate_arrays=args.store_rate_arrays,
            )
            full_cache[condition.case_id] = rows
            full_rows_all.extend(rows)
        species_scores = _score_species_by_trace(full_rows_all, list(full.species_names))
        score_source = "integrated_abs_rop_proxy"
        score_path = next((path for path in _score_candidates(assets_dir, benchmark) if path.exists()), None)
        if score_path is not None:
            reaction_scores = _read_score_csv(score_path, full.n_reactions)
            score_source = _rel(score_path)
        else:
            reaction_scores = _score_reactions_by_rop(
                full_mechanism,
                conditions,
                n_points=max(8, min(args.n_points, args.score_points)),
                surface_phase=surface_phase,
            )
        for method in methods:
            for level, ratio in effect_points:
                candidate_dir = out_dir / benchmark / method / level
                candidate_dir.mkdir(parents=True, exist_ok=True)
                try:
                    if method == "baseline":
                        candidate = _baseline_candidate(
                            benchmark,
                            full_mechanism,
                            candidate_dir,
                            level,
                            ratio,
                            conditions,
                            species_scores,
                            protected_species | set(major_species),
                            surface_phase,
                        )
                    elif method == "learnck":
                        candidate = _learnck_candidate(
                            benchmark,
                            full_mechanism,
                            candidate_dir,
                            level,
                            ratio,
                            reaction_scores,
                            score_source,
                            surface_phase,
                        )
                    elif method == "pooling_proxy":
                        candidate = _pooling_candidate(
                            benchmark,
                            full_mechanism,
                            candidate_dir,
                            level,
                            ratio,
                            conditions,
                            species_scores,
                            protected_species | set(major_species),
                            surface_phase,
                        )
                    else:  # pragma: no cover - guarded by parser
                        raise InputError(f"unknown method: {method}")
                    all_rows.extend(
                        _run_candidate(
                            benchmark,
                            full_mechanism,
                            conditions,
                            candidate,
                            candidate_dir,
                            full_cache=full_cache,
                            major_species=major_species,
                            n_points=args.n_points,
                            store_rate_arrays=args.store_rate_arrays,
                            max_network_reactions=args.max_network_reactions,
                            surface_phase=surface_phase,
                        )
                    )
                except Exception as exc:
                    status = {
                        "status": "failed",
                        "benchmark": benchmark,
                        "method": method,
                        "level": level,
                        "error": str(exc),
                    }
                    _write_status(candidate_dir / "status.json", status)
                    all_rows.append(
                        {
                            "benchmark": benchmark,
                            "method": method,
                            "method_label": method,
                            "level": level,
                            "case_id": "",
                            "status": "failed",
                            "error": str(exc),
                            "directory": _rel(candidate_dir),
                        }
                    )
    _write_csv(out_dir / "summary_all.csv", all_rows)
    _plot_effect_matrix(out_dir, all_rows)
    _plot_effect_curves(out_dir, all_rows)
    _write_index(out_dir, all_rows)
    print(f"wrote eval53 large reduction effects to {_rel(out_dir)}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", default="eval53_large")
    parser.add_argument("--methods", default="learnck,pooling_proxy")
    parser.add_argument("--levels", default="mild,medium,strong")
    parser.add_argument(
        "--sweep-ratios",
        default=None,
        help=(
            "Comma-separated keep ratios such as 1.0,0.9,0.75. "
            "When set, this overrides --levels and writes ratio_XXX directories."
        ),
    )
    parser.add_argument("--benchmarks", default=",".join(DEFAULT_BENCHMARKS))
    parser.add_argument(
        "--case-ids",
        default=None,
        help=(
            "Optional case filter. Use case_id for all benchmarks or "
            "benchmark:case_id for benchmark-specific cases."
        ),
    )
    parser.add_argument("--assets-dir", default="benchmarks/assets/eval53_large")
    parser.add_argument("--out", default="reports/eval53_large_reduction_effects")
    parser.add_argument("--n-points", type=int, default=120)
    parser.add_argument("--score-points", type=int, default=60)
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--major-species", default="auto")
    parser.add_argument("--protected-species", default="CO,CO2,H2O,H2,O2,N2")
    parser.add_argument("--max-network-reactions", type=int, default=40)
    parser.add_argument(
        "--store-rate-arrays",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Store rop_net_json and wdot_json arrays in trajectory CSV files.",
    )
    return parser


if __name__ == "__main__":
    sys.exit(run(build_parser().parse_args()))
