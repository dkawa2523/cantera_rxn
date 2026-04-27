from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "reports" / "eval53_method_compression_sweep_ac_sif4"
OUT = SWEEP / "report_assets"

METHOD_ORDER = ("baseline_proxy", "learnck_style_proxy", "pooling_proxy")
METHOD_DIR = {
    "baseline_proxy": "baseline",
    "learnck_style_proxy": "learnck",
    "pooling_proxy": "pooling_proxy",
}
METHOD_LABEL = {
    "baseline_proxy": "baseline",
    "learnck_style_proxy": "learnCK",
    "pooling_proxy": "pooling",
}
COLORS = {
    "baseline_proxy": "#4E79A7",
    "learnck_style_proxy": "#F28E2B",
    "pooling_proxy": "#59A14F",
}
SELECTED_RATIOS = ("ratio_1000", "ratio_750", "ratio_500", "ratio_200", "ratio_050")
RATIO_ORDER = (
    "ratio_1000",
    "ratio_900",
    "ratio_750",
    "ratio_600",
    "ratio_500",
    "ratio_400",
    "ratio_300",
    "ratio_200",
    "ratio_100",
    "ratio_050",
)
CASE_IDS = {"ac": "dilute_N2", "sif4": "add_O2"}
RATIO_TICK_LABEL = {
    "ratio_1000": "keep 1.00\n(no reduction)",
    "ratio_900": "keep 0.90",
    "ratio_750": "keep 0.75",
    "ratio_600": "keep 0.60",
    "ratio_500": "keep 0.50",
    "ratio_400": "keep 0.40",
    "ratio_300": "keep 0.30",
    "ratio_200": "keep 0.20",
    "ratio_100": "keep 0.10",
    "ratio_050": "keep 0.05\n(strongest)",
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _safe_float(row: dict[str, str], key: str, default: float = math.nan) -> float:
    try:
        value = row.get(key, "")
        if value == "":
            return default
        return float(value)
    except ValueError:
        return default


def _safe_stem(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in value).strip("_")


def _ratio_tick_label(level: str) -> str:
    return RATIO_TICK_LABEL.get(level, level.replace("ratio_", "keep "))


def _success_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row.get("status") == "success"]


def _case_key(row: dict[str, str]) -> tuple[str, str]:
    return row["benchmark"], row["case_id"]


def _species_error(row: dict[str, str]) -> float:
    total = 0.0
    for key, value in row.items():
        if key.startswith("X_") and key.endswith("_final_abs_delta") and value:
            total += abs(float(value))
    return total


def _save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        fig.tight_layout()
    except RuntimeError:
        pass
    fig.savefig(OUT / f"{stem}.png", dpi=180)
    fig.savefig(OUT / f"{stem}.svg")
    plt.close(fig)


def _read_cropped_image(path: Path) -> np.ndarray:
    image = mpimg.imread(path)
    rgb = image[..., :3]
    mask = np.any(rgb < 0.98, axis=2)
    if not mask.any():
        return image
    y_coords, x_coords = np.where(mask)
    pad = 12
    y0 = max(int(y_coords.min()) - pad, 0)
    y1 = min(int(y_coords.max()) + pad + 1, image.shape[0])
    x0 = max(int(x_coords.min()) - pad, 0)
    x1 = min(int(x_coords.max()) + pad + 1, image.shape[1])
    return image[y0:y1, x0:x1]


def plot_method_comparison(rows: list[dict[str, str]]) -> None:
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        subset = [
            row
            for row in _success_rows(rows)
            if row["benchmark"] == benchmark and row["case_id"] == case_id
        ]
        fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
        metrics = [
            ("T_final_rel_pct", "|final T relative error| [%]"),
            ("density_final_rel_pct", "|final density relative error| [%]"),
            ("T_max_delta", "|max T delta| [K]"),
            ("species_error", "sum |major species final error|"),
        ]
        for method in METHOD_ORDER:
            method_rows = [
                row for row in subset if row.get("method_label") == method
            ]
            method_rows.sort(key=lambda row: _safe_float(row, "reaction_reduction_pct"))
            if not method_rows:
                continue
            x = np.asarray([_safe_float(row, "reaction_reduction_pct") for row in method_rows])
            y_values = {
                "T_final_rel_pct": [abs(_safe_float(row, "T_final_rel_pct")) for row in method_rows],
                "density_final_rel_pct": [
                    abs(_safe_float(row, "density_final_rel_pct")) for row in method_rows
                ],
                "T_max_delta": [abs(_safe_float(row, "T_max_delta")) for row in method_rows],
                "species_error": [_species_error(row) for row in method_rows],
            }
            for ax, (metric, ylabel) in zip(axes.ravel(), metrics):
                ax.plot(
                    x,
                    y_values[metric],
                    marker="o",
                    linewidth=1.8,
                    label=METHOD_LABEL[method],
                    color=COLORS[method],
                )
                for row, xi, yi in zip(method_rows, x, y_values[metric]):
                    if row["level"] in SELECTED_RATIOS:
                        ax.annotate(
                            row["level"].replace("ratio_", ""),
                            (xi, yi),
                            xytext=(3, 4),
                            textcoords="offset points",
                            fontsize=7,
                        )
                ax.set_xlabel("actual reaction reduction [%]")
                ax.set_ylabel(ylabel)
                ax.grid(alpha=0.25)
        for ax in axes.ravel():
            ax.legend(frameon=False, fontsize=8)
        fig.suptitle(f"{benchmark} / {case_id}: method sweep Cantera error vs compression")
        _save(fig, f"method_comparison_{_safe_stem(benchmark)}_{_safe_stem(case_id)}")


def _trajectory_path(row: dict[str, str], kind: str) -> Path:
    return ROOT / row["directory"] / f"trajectories_{kind}.csv"


def _rows_by_method_level(rows: list[dict[str, str]], benchmark: str, case_id: str) -> dict[tuple[str, str], dict[str, str]]:
    return {
        (row.get("method_label", row.get("method", "")), row["level"]): row
        for row in _success_rows(rows)
        if row["benchmark"] == benchmark and row["case_id"] == case_id
    }


def _match_row(
    rows: list[dict[str, str]],
    benchmark: str,
    method: str,
    level: str,
    *,
    success_only: bool = False,
) -> dict[str, str] | None:
    raw_method = METHOD_DIR[method]
    for row in rows:
        if row.get("benchmark") != benchmark or row.get("level") != level:
            continue
        if success_only and row.get("status") != "success":
            continue
        label = row.get("method_label") or row.get("method")
        if label == method or row.get("method") == raw_method:
            return row
    return None


def _max_final_error(row: dict[str, str]) -> float:
    return max(
        abs(_safe_float(row, "T_final_rel_pct", 0.0)),
        abs(_safe_float(row, "density_final_rel_pct", 0.0)),
    )


def _is_accepted(row: dict[str, str]) -> bool:
    return (
        row.get("status") == "success"
        and abs(_safe_float(row, "T_final_rel_pct")) <= 1.0
        and abs(_safe_float(row, "density_final_rel_pct")) <= 1.0
    )


def _quality_color(error: float | None, failed: bool = False) -> str:
    if failed:
        return "#adb5bd"
    if error is None or math.isnan(error):
        return "#f1f3f5"
    if error <= 1.0:
        return "#2f9e44"
    if error <= 5.0:
        return "#ffd43b"
    if error <= 20.0:
        return "#f76707"
    return "#c92a2a"


def plot_quality_grids(rows: list[dict[str, str]]) -> None:
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        fig, ax = plt.subplots(figsize=(14.5, 4.8))
        ax.set_title(
            f"{benchmark} / {case_id}: target keep-ratio quality map "
            "(left=no reduction, right=stronger compression)"
        )
        ax.set_xlim(-0.5, len(RATIO_ORDER) - 0.5)
        ax.set_ylim(-0.5, len(METHOD_ORDER) - 0.5)
        ax.set_xticks(range(len(RATIO_ORDER)))
        ax.set_xticklabels([_ratio_tick_label(level) for level in RATIO_ORDER], rotation=35, ha="right")
        ax.set_yticks(range(len(METHOD_ORDER)))
        ax.set_yticklabels([METHOD_LABEL[method] for method in METHOD_ORDER])
        for y, method in enumerate(METHOD_ORDER):
            for x, level in enumerate(RATIO_ORDER):
                row = _match_row(rows, benchmark, method, level)
                if row is None:
                    color = _quality_color(None)
                    text = "n/a"
                elif row.get("status") != "success":
                    color = _quality_color(None, failed=True)
                    text = "failed"
                else:
                    error = _max_final_error(row)
                    color = _quality_color(error)
                    text = (
                        f"err {error:.1f}%\n"
                        f"Rred {_safe_float(row, 'reaction_reduction_pct'):.1f}%\n"
                        f"R {int(_safe_float(row, 'reactions_after', 0.0))}"
                    )
                rect = plt.Rectangle((x - 0.5, y - 0.5), 1.0, 1.0, facecolor=color, edgecolor="white")
                ax.add_patch(rect)
                text_color = "white" if color in {"#c92a2a", "#2f9e44"} else "#1f2933"
                ax.text(x, y, text, ha="center", va="center", fontsize=8, color=text_color)
        ax.set_xlabel("target keep ratio: fraction intended to remain; stronger compression to the right")
        ax.set_ylabel("method")
        ax.invert_yaxis()
        legend_items = [
            ("<=1% stable", "#2f9e44"),
            ("1-5% watch", "#ffd43b"),
            ("5-20% poor", "#f76707"),
            (">20% reject", "#c92a2a"),
            ("failed", "#adb5bd"),
        ]
        for index, (label, color) in enumerate(legend_items):
            ax.scatter([], [], marker="s", s=90, color=color, label=label)
        ax.legend(ncol=len(legend_items), frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2))
        _save(fig, f"quality_grid_{_safe_stem(benchmark)}_{_safe_stem(case_id)}")


def plot_cantera_delta_montage(rows: list[dict[str, str]]) -> None:
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        lookup = _rows_by_method_level(rows, benchmark, case_id)
        fig, axes = plt.subplots(len(METHOD_ORDER), 2, figsize=(12.0, 9.2), sharex=False)
        fig.suptitle(f"{benchmark} / {case_id}: reduced - full Cantera trajectory differences")
        for row_index, method in enumerate(METHOD_ORDER):
            ax_t = axes[row_index][0]
            ax_rho = axes[row_index][1]
            for level in SELECTED_RATIOS:
                row = lookup.get((method, level))
                if row is None:
                    continue
                full_path = _trajectory_path(row, "full")
                reduced_path = _trajectory_path(row, "reduced")
                if not full_path.exists() or not reduced_path.exists():
                    continue
                full = _read_csv(full_path)
                reduced = _read_csv(reduced_path)
                time = np.asarray([_safe_float(item, "time_s") for item in full])
                full_t = np.asarray([_safe_float(item, "T") for item in full])
                red_t = np.asarray([_safe_float(item, "T") for item in reduced])
                full_rho = np.asarray([_safe_float(item, "density") for item in full])
                red_rho = np.asarray([_safe_float(item, "density") for item in reduced])
                ax_t.plot(time, red_t - full_t, marker=None, linewidth=1.5, label=level.replace("ratio_", ""))
                ax_rho.plot(
                    time,
                    (red_rho - full_rho) / np.maximum(np.abs(full_rho), 1.0e-300) * 100.0,
                    marker=None,
                    linewidth=1.5,
                    label=level.replace("ratio_", ""),
                )
            ax_t.axhline(0.0, color="#202124", linewidth=0.8)
            ax_rho.axhline(0.0, color="#202124", linewidth=0.8)
            ax_t.set_title(f"{METHOD_LABEL[method]}: T error")
            ax_rho.set_title(f"{METHOD_LABEL[method]}: density error")
            ax_t.set_ylabel("T_red - T_full [K]")
            ax_rho.set_ylabel("density rel error [%]")
            ax_t.set_xlabel("time [s]")
            ax_rho.set_xlabel("time [s]")
            ax_t.grid(alpha=0.25)
            ax_rho.grid(alpha=0.25)
            ax_t.legend(frameon=False, fontsize=7, ncol=3)
            ax_rho.legend(frameon=False, fontsize=7, ncol=3)
        _save(fig, f"cantera_delta_montage_{_safe_stem(benchmark)}_{_safe_stem(case_id)}")


def plot_network_atlas(rows: list[dict[str, str]]) -> None:
    status_lookup = {
        (row["benchmark"], row.get("method_label", row.get("method", "")), row["level"]): row
        for row in rows
    }
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        fig, axes = plt.subplots(
            len(METHOD_ORDER),
            len(SELECTED_RATIOS),
            figsize=(17.0, 9.2),
        )
        fig.suptitle(f"{benchmark} / {case_id}: network after reduction by method and ratio")
        for row_index, method in enumerate(METHOD_ORDER):
            for col_index, level in enumerate(SELECTED_RATIOS):
                ax = axes[row_index][col_index]
                row = status_lookup.get((benchmark, method, level))
                if row is None:
                    raw = status_lookup.get((benchmark, method.replace("_proxy", ""), level))
                    row = raw
                title = f"{METHOD_LABEL[method]}\n{level.replace('ratio_', '')}"
                if row is None or row.get("status") != "success":
                    ax.set_facecolor("#f1f3f5")
                    ax.text(0.5, 0.5, "failed / no candidate", ha="center", va="center", fontsize=10)
                    ax.set_title(title)
                    ax.set_axis_off()
                    continue
                image_path = ROOT / row["directory"] / "network_after.png"
                if image_path.exists():
                    ax.imshow(_read_cropped_image(image_path))
                    title += f"\nR {row.get('reactions_after', '?')}"
                else:
                    ax.set_facecolor("#f1f3f5")
                    ax.text(0.5, 0.5, "no image", ha="center", va="center", fontsize=10)
                ax.set_title(title, fontsize=9)
                ax.set_axis_off()
        _save(fig, f"network_atlas_{_safe_stem(benchmark)}_{_safe_stem(case_id)}")


def _summary_row(
    summary: list[dict[str, object]],
    benchmark: str,
    case_id: str,
    method: str,
) -> dict[str, object] | None:
    method_name = METHOD_LABEL[method]
    for row in summary:
        if (
            row.get("benchmark") == benchmark
            and row.get("case_id") == case_id
            and row.get("method") == method_name
        ):
            return row
    return None


def _summary_level(
    summary: list[dict[str, object]],
    benchmark: str,
    case_id: str,
    method: str,
    key: str,
) -> str:
    row = _summary_row(summary, benchmark, case_id, method)
    if row is None:
        return ""
    value = row.get(key, "")
    return str(value) if value else ""


def plot_compression_delta_status(
    rows: list[dict[str, str]],
    summary: list[dict[str, object]],
) -> None:
    metrics = [
        ("T_final_rel_pct", "final T error [%]"),
        ("density_final_rel_pct", "final density error [%]"),
    ]
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        fig, axes = plt.subplots(len(METHOD_ORDER), len(metrics), figsize=(15.0, 9.8), sharex=True)
        fig.suptitle(
            (
                f"{benchmark} / {case_id}: reduced - full final differences by target keep ratio\n"
                "green circle=accept, yellow star=accepted boundary, blue circle=success only, gray x=failed"
            ),
            fontsize=13,
        )

        metric_limits: dict[str, float] = {}
        for metric, _label in metrics:
            values = [
                abs(_safe_float(row, metric))
                for row in rows
                if row.get("benchmark") == benchmark and row.get("status") == "success"
            ]
            metric_limits[metric] = max(2.0, max(values, default=2.0) * 1.18)

        for row_index, method in enumerate(METHOD_ORDER):
            boundary = _summary_level(summary, benchmark, case_id, method, "best_ratio_within_1pct")
            for col_index, (metric, ylabel) in enumerate(metrics):
                ax = axes[row_index][col_index]
                ax.axhspan(-1.0, 1.0, color="#d3f9d8", alpha=0.55, zorder=0)
                ax.axhline(0.0, color="#202124", linewidth=0.9, zorder=1)
                ax.axhline(1.0, color="#2f9e44", linewidth=0.8, linestyle=":", zorder=1)
                ax.axhline(-1.0, color="#2f9e44", linewidth=0.8, linestyle=":", zorder=1)
                success_x: list[int] = []
                success_y: list[float] = []
                failed_x: list[int] = []

                for x, level in enumerate(RATIO_ORDER):
                    row = _match_row(rows, benchmark, method, level)
                    if row is None or row.get("status") != "success":
                        failed_x.append(x)
                        continue
                    success_x.append(x)
                    success_y.append(_safe_float(row, metric))

                if success_x:
                    ax.plot(success_x, success_y, color="#778899", linewidth=1.3, alpha=0.55, zorder=2)

                y_limit = metric_limits[metric]
                failed_y = -y_limit * 0.91
                for x in failed_x:
                    ax.scatter(x, failed_y, marker="x", s=80, color="#868e96", linewidth=2.0, zorder=3)
                    ax.text(x, failed_y, "failed", ha="center", va="top", fontsize=7, color="#495057")

                for x, level in enumerate(RATIO_ORDER):
                    row = _match_row(rows, benchmark, method, level)
                    if row is None or row.get("status") != "success":
                        continue
                    y = _safe_float(row, metric)
                    if _is_accepted(row):
                        face = "#2f9e44"
                        edge = "#1b5e20"
                        accepted = True
                    else:
                        face = "#4c78a8"
                        edge = "#1f2933"
                        accepted = False
                    ax.scatter(
                        x,
                        y,
                        marker="o",
                        s=64,
                        facecolor=face,
                        edgecolor=edge,
                        linewidth=1.0,
                        zorder=4,
                    )
                    if level == boundary:
                        ax.scatter(
                            x,
                            y,
                            marker="*",
                            s=220,
                            facecolor="#ffd43b",
                            edgecolor="#202124",
                            linewidth=1.0,
                            zorder=5,
                        )
                        ax.annotate(
                            "ACCEPT\nboundary",
                            (x, y),
                            xytext=(6, 10),
                            textcoords="offset points",
                            fontsize=8,
                            color="#202124",
                            ha="left",
                            va="bottom",
                        )
                    elif not accepted and level in {"ratio_050", "ratio_200", "ratio_300"}:
                        ax.annotate(
                            level.replace("ratio_", ""),
                            (x, y),
                            xytext=(3, 5),
                            textcoords="offset points",
                            fontsize=7,
                            color="#1f2933",
                        )

                ax.set_ylim(-y_limit, y_limit)
                ax.set_xlim(-0.5, len(RATIO_ORDER) - 0.5)
                ax.set_ylabel(ylabel)
                ax.set_title(f"{METHOD_LABEL[method]}: {ylabel}", fontsize=10)
                ax.grid(axis="y", alpha=0.25)
                ax.set_xticks(range(len(RATIO_ORDER)))
                ax.set_xticklabels([_ratio_tick_label(level) for level in RATIO_ORDER], rotation=35, ha="right")
                if row_index == len(METHOD_ORDER) - 1:
                    ax.set_xlabel(
                        "target keep ratio: 1.00 keeps all reactions; 0.05 is strongest target compression"
                    )
        _save(fig, f"compression_delta_status_{_safe_stem(benchmark)}_{_safe_stem(case_id)}")


def plot_representative_trajectories(
    rows: list[dict[str, str]],
    summary: list[dict[str, object]],
) -> None:
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        fig, axes = plt.subplots(len(METHOD_ORDER), 2, figsize=(13.5, 10.0), sharex=False)
        fig.suptitle(
            f"{benchmark} / {case_id}: representative full vs reduced Cantera T trajectories",
            fontsize=14,
        )
        for row_index, method in enumerate(METHOD_ORDER):
            levels = [
                (
                    "accepted boundary",
                    _summary_level(summary, benchmark, case_id, method, "best_ratio_within_1pct"),
                ),
                (
                    "strongest success",
                    _summary_level(summary, benchmark, case_id, method, "max_success_ratio"),
                ),
            ]
            for col_index, (column_label, level) in enumerate(levels):
                ax = axes[row_index][col_index]
                ax.set_title(f"{METHOD_LABEL[method]}: {column_label}")
                if not level:
                    ax.set_facecolor("#f1f3f5")
                    ax.text(0.5, 0.5, "no accepted candidate", ha="center", va="center")
                    ax.set_axis_off()
                    continue
                row = _match_row(rows, benchmark, method, level, success_only=True)
                if row is None:
                    ax.set_facecolor("#f1f3f5")
                    ax.text(0.5, 0.5, f"{level}\nfailed / no data", ha="center", va="center")
                    ax.set_axis_off()
                    continue
                full_path = _trajectory_path(row, "full")
                reduced_path = _trajectory_path(row, "reduced")
                if not full_path.exists() or not reduced_path.exists():
                    ax.set_facecolor("#f1f3f5")
                    ax.text(0.5, 0.5, f"{level}\nmissing trajectory", ha="center", va="center")
                    ax.set_axis_off()
                    continue
                full = _read_csv(full_path)
                reduced = _read_csv(reduced_path)
                time = np.asarray([_safe_float(item, "time_s") for item in full])
                full_t = np.asarray([_safe_float(item, "T") for item in full])
                red_t = np.asarray([_safe_float(item, "T") for item in reduced])
                ax.plot(time, full_t, color="#202124", linewidth=2.2, label="full")
                ax.plot(
                    time,
                    red_t,
                    color=COLORS[method],
                    linewidth=2.0,
                    linestyle="--",
                    label="reduced",
                )
                ax.grid(alpha=0.25)
                ax.set_xlabel("time [s]")
                ax.set_ylabel("T [K]")
                ax.legend(frameon=False, fontsize=8, loc="best")
                text = (
                    f"{level.replace('ratio_', 'ratio ')}\n"
                    f"reaction reduction {_safe_float(row, 'reaction_reduction_pct'):.1f}%\n"
                    f"final T err {_safe_float(row, 'T_final_rel_pct'):.2f}%\n"
                    f"final density err {_safe_float(row, 'density_final_rel_pct'):.2f}%"
                )
                ax.text(
                    0.03,
                    0.97,
                    text,
                    transform=ax.transAxes,
                    va="top",
                    ha="left",
                    fontsize=8,
                    bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "#d9e2ec"},
                )
        _save(fig, f"trajectory_representatives_{_safe_stem(benchmark)}_{_safe_stem(case_id)}")


def plot_network_focus(
    rows: list[dict[str, str]],
    summary: list[dict[str, object]],
) -> None:
    columns = [
        ("no reduction", "ratio_1000"),
        ("accepted boundary", "best_ratio_within_1pct"),
        ("strongest success", "max_success_ratio"),
    ]
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        fig, axes = plt.subplots(len(METHOD_ORDER), len(columns), figsize=(14.5, 11.0))
        fig.suptitle(
            f"{benchmark} / {case_id}: focused network states for compression decisions",
            fontsize=14,
        )
        for row_index, method in enumerate(METHOD_ORDER):
            for col_index, (column_label, value) in enumerate(columns):
                ax = axes[row_index][col_index]
                level = value
                if value != "ratio_1000":
                    level = _summary_level(summary, benchmark, case_id, method, value)
                title = f"{METHOD_LABEL[method]}\n{column_label}"
                if not level:
                    ax.set_facecolor("#f1f3f5")
                    ax.text(0.5, 0.5, "no accepted candidate", ha="center", va="center", fontsize=11)
                    ax.set_title(title, fontsize=10)
                    ax.set_axis_off()
                    continue
                row = _match_row(rows, benchmark, method, level, success_only=True)
                if row is None:
                    ax.set_facecolor("#f1f3f5")
                    ax.text(0.5, 0.5, f"{level}\nfailed / no data", ha="center", va="center", fontsize=11)
                    ax.set_title(title, fontsize=10)
                    ax.set_axis_off()
                    continue
                image_path = ROOT / row["directory"] / "network_after.png"
                reactions_after = int(round(_safe_float(row, "reactions_after", 0.0)))
                if reactions_after <= 0:
                    ax.set_facecolor("#f8f9fa")
                    ax.text(0.5, 0.5, "0 reactions\nnetwork collapsed", ha="center", va="center", fontsize=11)
                elif image_path.exists():
                    ax.imshow(_read_cropped_image(image_path))
                else:
                    ax.set_facecolor("#f1f3f5")
                    ax.text(0.5, 0.5, "no image", ha="center", va="center", fontsize=11)
                ax.set_title(
                    (
                        f"{title}\n{level.replace('ratio_', 'ratio ')}"
                        f" | R {row.get('reactions_after', '?')}"
                        f" | Rred {_safe_float(row, 'reaction_reduction_pct'):.1f}%"
                    ),
                    fontsize=9,
                )
                ax.set_axis_off()
        _save(fig, f"network_focus_{_safe_stem(benchmark)}_{_safe_stem(case_id)}")


def plot_network_decision_focus(
    rows: list[dict[str, str]],
    summary: list[dict[str, object]],
) -> None:
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        learnck_accept = _summary_level(
            summary, benchmark, case_id, "learnck_style_proxy", "best_ratio_within_1pct"
        )
        learnck_strong = _summary_level(
            summary, benchmark, case_id, "learnck_style_proxy", "max_success_ratio"
        )
        pooling_strong = _summary_level(
            summary, benchmark, case_id, "pooling_proxy", "max_success_ratio"
        )
        panels = [
            ("reference network", "baseline_proxy", "ratio_1000"),
            ("accepted learnCK network", "learnck_style_proxy", learnck_accept),
            ("over-compressed learnCK network", "learnck_style_proxy", learnck_strong),
            ("pooling degenerate network", "pooling_proxy", pooling_strong),
        ]
        fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.2))
        fig.suptitle(
            f"{benchmark} / {case_id}: four network states that explain the decision",
            fontsize=14,
        )
        for ax, (panel_label, method, level) in zip(axes.ravel(), panels):
            title = f"{panel_label}\n{METHOD_LABEL[method]}"
            if not level:
                ax.set_facecolor("#f1f3f5")
                ax.text(0.5, 0.5, "no candidate", ha="center", va="center", fontsize=12)
                ax.set_title(title, fontsize=11)
                ax.set_axis_off()
                continue
            row = _match_row(rows, benchmark, method, level, success_only=True)
            if row is None:
                ax.set_facecolor("#f1f3f5")
                ax.text(0.5, 0.5, f"{level}\nfailed / no data", ha="center", va="center", fontsize=12)
                ax.set_title(title, fontsize=11)
                ax.set_axis_off()
                continue
            image_path = ROOT / row["directory"] / "network_after.png"
            reactions_after = int(round(_safe_float(row, "reactions_after", 0.0)))
            if reactions_after <= 0:
                ax.set_facecolor("#f8f9fa")
                ax.text(0.5, 0.5, "0 reactions\nnetwork collapsed", ha="center", va="center", fontsize=13)
            elif image_path.exists():
                ax.imshow(_read_cropped_image(image_path))
            else:
                ax.set_facecolor("#f1f3f5")
                ax.text(0.5, 0.5, "no image", ha="center", va="center", fontsize=12)
            ax.set_title(
                (
                    f"{title}\n{level.replace('ratio_', 'ratio ')}"
                    f" | reactions {row.get('reactions_after', '?')}"
                    f" | reduction {_safe_float(row, 'reaction_reduction_pct'):.1f}%"
                    f" | err {_max_final_error(row):.1f}%"
                ),
                fontsize=10,
            )
            ax.set_axis_off()
        _save(fig, f"network_decision_focus_{_safe_stem(benchmark)}_{_safe_stem(case_id)}")


def write_decision_summary(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for benchmark, case_id in sorted({_case_key(row) for row in _success_rows(rows)}):
        for method in METHOD_ORDER:
            successes = [
                row
                for row in _success_rows(rows)
                if row["benchmark"] == benchmark
                and row["case_id"] == case_id
                and row.get("method_label") == method
            ]
            successes.sort(key=lambda row: _safe_float(row, "reaction_reduction_pct"))
            failed = [
                row
                for row in rows
                if row["benchmark"] == benchmark
                and row.get("method_label", row.get("method")) in {method, method.replace("_proxy", "")}
                and row.get("status") != "success"
            ]
            within_1pct = [
                row
                for row in successes
                if abs(_safe_float(row, "T_final_rel_pct")) <= 1.0
                and abs(_safe_float(row, "density_final_rel_pct")) <= 1.0
            ]
            within_5pct = [
                row
                for row in successes
                if abs(_safe_float(row, "T_final_rel_pct")) <= 5.0
                and abs(_safe_float(row, "density_final_rel_pct")) <= 5.0
            ]
            best_1pct = max(within_1pct, key=lambda row: _safe_float(row, "reaction_reduction_pct"), default=None)
            best_5pct = max(within_5pct, key=lambda row: _safe_float(row, "reaction_reduction_pct"), default=None)
            max_success = max(successes, key=lambda row: _safe_float(row, "reaction_reduction_pct"), default=None)
            output.append(
                {
                    "benchmark": benchmark,
                    "case_id": case_id,
                    "method": METHOD_LABEL[method],
                    "success_count": len(successes),
                    "failed_count": len(failed),
                    "best_ratio_within_1pct": best_1pct["level"] if best_1pct else "",
                    "reaction_reduction_within_1pct": _safe_float(best_1pct, "reaction_reduction_pct", "") if best_1pct else "",
                    "best_ratio_within_5pct": best_5pct["level"] if best_5pct else "",
                    "reaction_reduction_within_5pct": _safe_float(best_5pct, "reaction_reduction_pct", "") if best_5pct else "",
                    "max_success_ratio": max_success["level"] if max_success else "",
                    "max_success_reaction_reduction": _safe_float(max_success, "reaction_reduction_pct", "") if max_success else "",
                    "max_success_T_final_rel_pct": _safe_float(max_success, "T_final_rel_pct", "") if max_success else "",
                    "max_success_density_final_rel_pct": _safe_float(max_success, "density_final_rel_pct", "") if max_success else "",
                }
            )
    _write_csv(OUT / "method_sweep_decision_summary.csv", output)
    return output


def plot_decision_dashboard(summary: list[dict[str, object]]) -> None:
    cases = [("ac", "dilute_N2"), ("sif4", "add_O2")]
    methods = ("baseline", "learnCK", "pooling")
    lookup = {
        (str(row["benchmark"]), str(row["case_id"]), str(row["method"])): row
        for row in summary
    }
    usable = np.full((len(methods), len(cases)), np.nan)
    max_error = np.full((len(methods), len(cases)), np.nan)
    labels: list[list[str]] = [["" for _case in cases] for _method in methods]
    for row_index, method in enumerate(methods):
        for col_index, (benchmark, case_id) in enumerate(cases):
            row = lookup.get((benchmark, case_id, method))
            if row is None:
                continue
            reduction = row.get("reaction_reduction_within_1pct", "")
            if reduction != "":
                usable[row_index, col_index] = float(reduction)
            t_err = abs(float(row.get("max_success_T_final_rel_pct") or 0.0))
            rho_err = abs(float(row.get("max_success_density_final_rel_pct") or 0.0))
            max_error[row_index, col_index] = max(t_err, rho_err)
            labels[row_index][col_index] = (
                f"{row.get('best_ratio_within_1pct') or 'none'}\n"
                f"fail {row.get('failed_count')}"
            )

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    fig.suptitle("Method compression sweep: decision dashboard", fontsize=14)

    im0 = axes[0].imshow(usable, cmap="YlGn", vmin=0.0, vmax=100.0)
    axes[0].set_title("Largest reaction reduction with <=1% final T and density error")
    axes[0].set_xticks(range(len(cases)))
    axes[0].set_xticklabels([f"{b}\n{c}" for b, c in cases])
    axes[0].set_yticks(range(len(methods)))
    axes[0].set_yticklabels(methods)
    for i in range(len(methods)):
        for j in range(len(cases)):
            value = usable[i, j]
            text = "none" if not np.isfinite(value) else f"{value:.1f}%"
            axes[0].text(j, i, f"{text}\n{labels[i][j]}", ha="center", va="center", fontsize=9)
    cbar0 = fig.colorbar(im0, ax=axes[0], shrink=0.85)
    cbar0.set_label("reaction reduction [%]")

    im1 = axes[1].imshow(max_error, cmap="RdYlGn_r", vmin=0.0, vmax=50.0)
    axes[1].set_title("Worst final error at strongest successful candidate")
    axes[1].set_xticks(range(len(cases)))
    axes[1].set_xticklabels([f"{b}\n{c}" for b, c in cases])
    axes[1].set_yticks(range(len(methods)))
    axes[1].set_yticklabels(methods)
    for i in range(len(methods)):
        for j in range(len(cases)):
            value = max_error[i, j]
            axes[1].text(j, i, f"{value:.1f}%", ha="center", va="center", fontsize=10)
    cbar1 = fig.colorbar(im1, ax=axes[1], shrink=0.85)
    cbar1.set_label("max(|T err|, |density err|) [%]")
    _save(fig, "method_sweep_decision_dashboard")


def write_index(summary: list[dict[str, object]]) -> None:
    lines = [
        "# Method Compression Sweep Report Assets",
        "",
        "- [method_sweep_decision_summary.csv](method_sweep_decision_summary.csv)",
        "",
        "## Figures",
        "",
    ]
    for path in sorted(OUT.glob("*.png")):
        lines.append(f"- [{path.name}]({path.name})")
    lines.extend(["", "## Decision Summary", ""])
    lines.append("| benchmark | case | method | success | failed | best <=1% | best <=5% |")
    lines.append("|---|---|---|---:|---:|---:|---:|")
    for row in summary:
        lines.append(
            f"| {row['benchmark']} | {row['case_id']} | {row['method']} | "
            f"{row['success_count']} | {row['failed_count']} | "
            f"{row['best_ratio_within_1pct']} | {row['best_ratio_within_5pct']} |"
        )
    (OUT / "index.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    rows = _read_csv(SWEEP / "summary_all.csv")
    OUT.mkdir(parents=True, exist_ok=True)
    plot_method_comparison(rows)
    plot_quality_grids(rows)
    plot_cantera_delta_montage(rows)
    plot_network_atlas(rows)
    summary = write_decision_summary(rows)
    plot_decision_dashboard(summary)
    plot_compression_delta_status(rows, summary)
    plot_representative_trajectories(rows, summary)
    plot_network_focus(rows, summary)
    plot_network_decision_focus(rows, summary)
    write_index(summary)
    print(f"wrote method sweep report assets to {OUT}")


if __name__ == "__main__":
    main()
