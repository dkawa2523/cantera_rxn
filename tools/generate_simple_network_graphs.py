from __future__ import annotations

import argparse
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


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

NODE_RE = re.compile(r"^\s*([A-Za-z0-9_]+)\s+\[(.*)\];\s*$")
EDGE_RE = re.compile(r"^\s*([A-Za-z0-9_]+)\s*->\s*([A-Za-z0-9_]+)\s+\[(.*)\];\s*$")
ATTR_RE_TEMPLATE = r'{name}="((?:\\.|[^"\\])*)"'

PHASE_COLORS = {
    "gas": "#4E79A7",
    "surface": "#F28E2B",
    "bulk": "#59A14F",
    "other": "#59A14F",
}


@dataclass(frozen=True)
class Node:
    node_id: str
    name: str
    phase: str
    members: tuple[str, ...]


@dataclass(frozen=True)
class Edge:
    src: str
    dst: str
    reaction_ids: str
    equation: str


def _attr(attrs: str, name: str) -> str:
    match = re.search(ATTR_RE_TEMPLATE.format(name=re.escape(name)), attrs)
    return match.group(1) if match else ""


def _unescape_label(value: str) -> str:
    value = value.replace(r"\"", '"')
    value = value.replace(r"\\", "\\")
    return value.split(r"\n", 1)[0].strip()


def _phase_from_attrs(attrs: str) -> str:
    fill = _attr(attrs, "fillcolor").upper()
    tooltip = _attr(attrs, "tooltip").lower()
    if fill == "#D9ECFF":
        return "gas"
    if fill == "#FCE1C2":
        return "surface"
    if fill == "#D9F2D9":
        return "bulk"
    if "phase: gas" in tooltip or "gas:" in tooltip:
        return "gas"
    if "phase: surface" in tooltip or "surface:" in tooltip:
        return "surface"
    if "phase: bulk" in tooltip or "bulk:" in tooltip:
        return "bulk"
    return "other"


def _tooltip_field(tooltip: str, name: str) -> str:
    prefix = f"{name}:"
    for part in tooltip.replace("\\n", "\n").splitlines():
        if part.startswith(prefix):
            return (
                part[len(prefix) :].strip()
                .replace("…", "...")
                .replace("窶ｦ", "...")
                .replace("乧", "...")
            )
    return ""


def _members_from_attrs(attrs: str, label: str) -> tuple[str, ...]:
    tooltip = _attr(attrs, "tooltip")
    raw = _tooltip_field(tooltip, "members")
    if not raw:
        raw = _tooltip_field(tooltip, "name")
    if not raw:
        return (label,)
    members = tuple(item.strip() for item in raw.split(",") if item.strip())
    return members or (label,)


def parse_dot(path: Path) -> tuple[list[Node], list[Edge]]:
    nodes: dict[str, Node] = {}
    edges: list[Edge] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        edge_match = EDGE_RE.match(line)
        if edge_match:
            src, dst, attrs = edge_match.groups()
            tooltip = _attr(attrs, "tooltip")
            reaction_ids = _tooltip_field(tooltip, "reactions")
            equation = _tooltip_field(tooltip, "equations")
            edges.append(Edge(src=src, dst=dst, reaction_ids=reaction_ids, equation=equation))
            continue
        node_match = NODE_RE.match(line)
        if not node_match:
            continue
        node_id, attrs = node_match.groups()
        if node_id in {"graph", "node", "edge"}:
            continue
        if node_id.startswith("lane_") or 'shape="plaintext"' in attrs:
            continue
        label = _unescape_label(_attr(attrs, "label")) or node_id
        nodes[node_id] = Node(
            node_id=node_id,
            name=label,
            phase=_phase_from_attrs(attrs),
            members=_members_from_attrs(attrs, label),
        )
    edges = [edge for edge in edges if edge.src in nodes and edge.dst in nodes]
    return list(nodes.values()), edges


def _edge_color(src_phase: str, dst_phase: str) -> str:
    if src_phase == "gas" and dst_phase == "gas":
        return PHASE_COLORS["gas"]
    if src_phase == "surface" and dst_phase == "surface":
        return PHASE_COLORS["surface"]
    return PHASE_COLORS["bulk"]


def write_simple_dot(
    path: Path,
    *,
    title: str,
    nodes: list[Node],
    edges: list[Edge],
    summary_species: int,
    summary_reactions: int,
    node_numbers: dict[str, str] | None = None,
) -> None:
    phase_by_id = {node.node_id: node.phase for node in nodes}
    node_numbers = node_numbers or {}
    lines = [
        "digraph G {",
        '  graph [layout="sfdp", overlap="prism", splines="true", outputorder="edgesfirst", bgcolor="white", pad="0.25", labelloc="t", labeljust="l"];',
        f'  label="{title}";',
        '  node [shape="circle", fixedsize="true", width="0.30", height="0.30", label="", style="filled", penwidth="0.6", color="#263238", fontname="Helvetica", fontsize="7"];',
        '  edge [arrowsize="0.35", penwidth="0.55", color="#77777788"];',
    ]
    for node in nodes:
        color = PHASE_COLORS.get(node.phase, PHASE_COLORS["other"])
        tooltip = "; ".join(node.members).replace("\\", "\\\\").replace('"', r"\"")
        label = node_numbers.get(node.node_id, "")
        lines.append(f'  {node.node_id} [fillcolor="{color}", label="{label}", tooltip="{tooltip}"];')
    for edge in edges:
        color = _edge_color(phase_by_id.get(edge.src, "other"), phase_by_id.get(edge.dst, "other"))
        lines.append(f'  {edge.src} -> {edge.dst} [color="{color}88"];')
    lines.append("}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _natural_node_key(node: Node) -> tuple[str, int]:
    match = re.match(r"([A-Za-z_]+)(\d+)$", node.node_id)
    if match:
        return (match.group(1), int(match.group(2)))
    return (node.node_id, 0)


def _shorten(text: str, limit: int = 180) -> str:
    text = " ".join(text.split())
    if len(text) <= limit:
        return text
    return text[: limit - 3].rstrip() + "..."


def _format_members(members: tuple[str, ...]) -> str:
    text = ", ".join(members)
    if len(members) > 1:
        text = f"{len(members)} merged: {text}"
    return text


def write_node_map(
    path: Path,
    *,
    title: str,
    nodes: list[Node],
    edges: list[Edge],
    summary_species: int,
    summary_reactions: int,
    node_numbers: dict[str, str],
    png_name: str,
    svg_name: str,
) -> None:
    incident: dict[str, list[str]] = {node.node_id: [] for node in nodes}
    for edge in edges:
        rid = edge.reaction_ids or "?"
        eq = edge.equation or f"{edge.src} -> {edge.dst}"
        item = f"R{rid}: {eq}"
        incident.setdefault(edge.src, []).append(item)
        incident.setdefault(edge.dst, []).append(item)

    lines = [
        f"# {title}",
        "",
        f"![graph]({png_name})",
        "",
        f"- summary nodes: {summary_species}",
        f"- summary reactions: {summary_reactions}",
        f"- drawn nodes: {len(nodes)}",
        f"- drawn edges: {len(edges)}",
        "- colors: gas=blue, surface=orange, bulk/mixed=green",
        "",
    ]
    for node in sorted(nodes, key=_natural_node_key):
        label = node_numbers[node.node_id]
        color = {"gas": "blue", "surface": "orange", "bulk": "green"}.get(node.phase, "green")
        lines.extend(
            [
                f"## {label} ({color})",
                "",
                f"Names: {_format_members(node.members)}",
                "",
                "Reactions:",
            ]
        )
        reactions = list(dict.fromkeys(incident.get(node.node_id, [])))
        if reactions:
            lines.extend(f"- {reaction}" for reaction in reactions)
        else:
            lines.append("- none")
        lines.append("")
    lines.extend(["", f"SVG: [{svg_name}]({svg_name})", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def render(dot_cmd: str, dot_path: Path) -> None:
    for fmt in ("svg", "png"):
        out_path = dot_path.with_suffix(f".{fmt}")
        subprocess.run(
            [dot_cmd, f"-T{fmt}", "-Ksfdp", str(dot_path), "-o", str(out_path)],
            check=True,
        )


def _entry_by_mode(summary_path: Path, mode: str) -> dict:
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    for entry in data["entries"]:
        if entry["mode"] == mode:
            return entry
    raise KeyError(f"mode not found in {summary_path}: {mode}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--out-dir", default="reports/eval53viz_simple_networks")
    parser.add_argument("--dot-cmd", default="dot")
    parser.add_argument("--with-node-map", action="store_true")
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for benchmark, mode, summary_rel in RUNS:
        entry = _entry_by_mode(root / summary_rel, mode)
        run_id = entry["run_id"]
        viz_dir = root / "reports" / run_id / "viz"
        variants = [
            (
                "before_raw",
                viz_dir / "network_initial_reaction_edge_detail.dot",
                int(entry["species_before"]),
                int(entry["reactions_before"]),
            ),
            (
                "after_reduced",
                viz_dir / "network_reaction_edge_detail.dot",
                int(entry["species_after"]),
                int(entry["reactions_after"]),
            ),
        ]
        for variant, source_dot, summary_species, summary_reactions in variants:
            nodes, edges = parse_dot(source_dot)
            node_numbers = {
                node.node_id: f"N{i + 1}"
                for i, node in enumerate(sorted(nodes, key=_natural_node_key))
            }
            stem = f"{run_id}_{variant}_simple"
            dot_path = out_dir / f"{stem}.dot"
            title = f"{benchmark} {mode} {variant.replace('_', ' ')}"
            write_simple_dot(
                dot_path,
                title=title,
                nodes=nodes,
                edges=edges,
                summary_species=summary_species,
                summary_reactions=summary_reactions,
                node_numbers=node_numbers if args.with_node_map else None,
            )
            render(args.dot_cmd, dot_path)
            map_name = ""
            if args.with_node_map:
                map_path = out_dir / f"{stem}_node_map.md"
                map_name = map_path.name
                write_node_map(
                    map_path,
                    title=title,
                    nodes=nodes,
                    edges=edges,
                    summary_species=summary_species,
                    summary_reactions=summary_reactions,
                    node_numbers=node_numbers,
                    png_name=dot_path.with_suffix(".png").name,
                    svg_name=dot_path.with_suffix(".svg").name,
                )
            rows.append(
                {
                    "benchmark": benchmark,
                    "mode": mode,
                    "variant": variant,
                    "summary_nodes": summary_species,
                    "summary_reactions": summary_reactions,
                    "drawn_nodes": len(nodes),
                    "drawn_edges": len(edges),
                    "svg": dot_path.with_suffix(".svg").name,
                    "png": dot_path.with_suffix(".png").name,
                    "map": map_name,
                }
            )

    index_lines = [
        "# Eval53viz Simple State-Reaction Networks",
        "",
        "Colors: gas blue, surface orange, bulk/mixed green. Node size is fixed; graph labels are short node IDs.",
        "",
        "| benchmark | mode | graph | summary nodes | summary reactions | drawn nodes | drawn edges | svg | png | node map |",
        "|---|---|---|---:|---:|---:|---:|---|---|---|",
    ]
    for row in rows:
        index_lines.append(
            "| {benchmark} | {mode} | {variant} | {summary_nodes} | {summary_reactions} | "
            "{drawn_nodes} | {drawn_edges} | [{svg}]({svg}) | [{png}]({png}) | {map_link} |".format(
                **row,
                map_link=f"[map]({row['map']})" if row.get("map") else "-",
            )
        )
    (out_dir / "index.md").write_text("\n".join(index_lines) + "\n", encoding="utf-8")
    print(f"wrote {len(rows)} graphs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
