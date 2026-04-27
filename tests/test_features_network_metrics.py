from __future__ import annotations

import json
import math
from typing import Any

import rxn_platform.backends.dummy  # noqa: F401
import rxn_platform.tasks.features  # noqa: F401
import rxn_platform.tasks.sim  # noqa: F401

from rxn_platform.core import ArtifactManifest
from rxn_platform.io_utils import write_json_atomic
from rxn_platform.registry import get
from rxn_platform.store import ArtifactStore


def _read_features(path) -> list[dict[str, Any]]:
    try:
        import pandas as pd
    except ImportError:
        pd = None
    if pd is not None:
        try:
            frame = pd.read_parquet(path)
            return frame.to_dict(orient="records")
        except Exception:
            pass
    try:
        import pyarrow.parquet as pq
    except ImportError:
        pq = None
    if pq is not None:
        try:
            table = pq.read_table(path)
            return table.to_pylist()
        except Exception:
            pass
    payload = json.loads(path.read_text(encoding="utf-8"))
    return list(payload.get("rows", []))


def _make_graph_manifest(graph_id: str) -> ArtifactManifest:
    return ArtifactManifest(
        schema_version=1,
        kind="graphs",
        id=graph_id,
        created_at="2026-01-18T00:00:00Z",
        parents=[],
        inputs={"graph": "demo"},
        config={"source": "test"},
        code={"version": "0.0.0"},
        provenance={"python": "3.11"},
    )


def _write_demo_graph(store: ArtifactStore, graph_id: str) -> None:
    payload = {
        "directed": True,
        "multigraph": False,
        "graph": {"name": "demo"},
        "nodes": [
            {"id": "A", "kind": "species"},
            {"id": "B", "kind": "species"},
            {"id": "C", "kind": "species"},
            {"id": "D", "kind": "species"},
        ],
        "links": [
            {"source": "A", "target": "B"},
            {"source": "A", "target": "C"},
            {"source": "B", "target": "C"},
            {"source": "C", "target": "D"},
        ],
    }
    manifest = _make_graph_manifest(graph_id)

    def _writer(base_dir):
        (base_dir / "graph.json").write_text(
            json.dumps(payload, ensure_ascii=True, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    store.ensure(manifest, writer=_writer)


def _write_gnn_dataset_with_extra_run(store: ArtifactStore, dataset_id: str) -> None:
    manifest = ArtifactManifest(
        schema_version=1,
        kind="gnn_datasets",
        id=dataset_id,
        created_at="2026-01-18T00:00:00Z",
        parents=[],
        inputs={"source": "test"},
        config={"source": "test"},
        code={"version": "0.0.0"},
        provenance={"python": "3.11"},
    )

    def _writer(base_dir):
        write_json_atomic(
            base_dir / "dataset.json",
            {
                "source": {"run_ids": ["known-run", "extra-run"], "graph_id": "g"},
                "windows": [
                    {"index": 0, "window": {"start_idx": 0, "end_idx": 1}},
                ],
                "nodes": {
                    "order": ["species_A"],
                    "meta": [
                        {"id": "species_A", "species": "A", "species_index": 0},
                    ],
                },
                "files": {"data_json": "data.json"},
            },
        )
        write_json_atomic(
            base_dir / "data.json",
            {
                "items": [
                    {
                        "run_id": "known-run",
                        "window_id": 0,
                        "x": [[2.0]],
                        "edge_index": [[], []],
                        "edge_attr": [],
                    },
                    {
                        "run_id": "extra-run",
                        "window_id": 0,
                        "x": [[9.0]],
                        "edge_index": [[], []],
                        "edge_attr": [],
                    },
                ]
            },
        )

    store.ensure(manifest, writer=_writer)


def _run_dummy_sim(store: ArtifactStore, temperature: float) -> str:
    sim_task = get("task", "sim.run")
    sim_cfg = {
        "sim": {
            "name": "dummy",
            "time": {"start": 0.0, "stop": 1.0, "steps": 2},
            "initial": {"T": temperature, "P": 1000.0},
            "ramp": {"T": 0.0, "P": 0.0},
            "species": ["A", "B"],
        }
    }
    result = sim_task(sim_cfg, store=store)
    return result.manifest.id


def test_network_metrics_features_from_graph(tmp_path) -> None:
    store = ArtifactStore(tmp_path / "artifacts")
    graph_id = "graph-demo"
    _write_demo_graph(store, graph_id)

    run_id = _run_dummy_sim(store, temperature=100.0)

    feat_task = get("task", "features.run")
    feat_cfg = {
        "inputs": {"run_id": run_id},
        "params": {
            "features": [
                {
                    "name": "network_metrics",
                    "params": {
                        "graph_id": graph_id,
                        "metrics": ["degree", "degree_centrality"],
                        "top_n": 2,
                    },
                }
            ]
        },
    }
    feat_result = feat_task(feat_cfg, store=store)

    rows = _read_features(feat_result.path / "features.parquet")
    by_name = {row["feature"]: row for row in rows}

    assert "network.degree.C" in by_name
    assert "network.degree.A" in by_name
    assert math.isclose(by_name["network.degree.C"]["value"], 3.0, rel_tol=0.0, abs_tol=1.0e-9)

    meta = json.loads(by_name["network.degree.C"]["meta_json"])
    assert meta.get("feature_kind") == "network_metrics"
    assert meta.get("metric") == "degree"
    assert meta.get("graph_id") == graph_id


def test_network_metrics_rank_stability_meta(tmp_path) -> None:
    store = ArtifactStore(tmp_path / "artifacts")
    graph_id = "graph-demo"
    _write_demo_graph(store, graph_id)

    run_id_a = _run_dummy_sim(store, temperature=100.0)
    run_id_b = _run_dummy_sim(store, temperature=110.0)

    feat_task = get("task", "features.run")
    feat_cfg = {
        "inputs": {"run_ids": [run_id_a, run_id_b]},
        "params": {
            "features": [
                {
                    "name": "network_metrics",
                    "params": {
                        "graph_id": graph_id,
                        "metrics": ["degree"],
                        "top_n": 3,
                    },
                }
            ]
        },
    }
    feat_result = feat_task(feat_cfg, store=store)

    rows = _read_features(feat_result.path / "features.parquet")
    target_row = next(
        row
        for row in rows
        if row["run_id"] == run_id_a and row["feature"] == "network.degree.C"
    )
    meta = json.loads(target_row["meta_json"])
    stability = meta.get("rank_stability")
    assert stability is not None
    assert stability.get("status") == "computed"
    assert math.isclose(stability.get("spearman_mean"), 1.0, rel_tol=0.0, abs_tol=1.0e-9)
    assert math.isclose(
        stability.get("top_k_jaccard_mean"), 1.0, rel_tol=0.0, abs_tol=1.0e-9
    )


def test_gnn_importance_respects_explicit_run_filter(tmp_path) -> None:
    store = ArtifactStore(tmp_path / "artifacts")
    dataset_id = "gnn-filter-demo"
    _write_gnn_dataset_with_extra_run(store, dataset_id)

    task = get("task", "features.gnn_importance")
    result = task(
        {
            "inputs": {
                "dataset_id": dataset_id,
                "graph_id": "g",
                "run_ids": ["known-run"],
            },
            "params": {
                "output": {
                    "include_species": True,
                    "include_reactions": False,
                }
            },
        },
        store=store,
    )

    rows = _read_features(result.path / "features.parquet")
    assert result.manifest.inputs["runs"] == ["known-run"]
    assert {row["run_id"] for row in rows} == {"known-run"}
    assert all(row["feature"] == "gnn_species_importance" for row in rows)
