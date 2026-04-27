import json

import rxn_platform.backends.dummy  # Register backend.
import rxn_platform.tasks.sim  # Register task.

from rxn_platform.backends.dummy import DummyBackend
from rxn_platform.registry import get, register, resolve_backend
from rxn_platform.store import ArtifactStore


def test_dummy_backend_generates_run_dataset() -> None:
    backend = DummyBackend()
    cfg = {
        "time": {"start": 0.0, "stop": 1.0, "steps": 3},
        "initial": {"T": 900.0, "P": 101325.0},
        "ramp": {"T": 5.0, "P": 100.0},
        "species": ["A", "B"],
    }
    dataset = backend.run(cfg)

    assert dataset.coords["time"]["data"] == [0.0, 0.5, 1.0]
    assert dataset.coords["species"]["data"] == ["A", "B"]
    assert dataset.data_vars["T"]["data"] == [900.0, 905.0, 910.0]
    assert dataset.data_vars["P"]["data"] == [101325.0, 101425.0, 101525.0]
    assert len(dataset.data_vars["X"]["data"]) == 3


def test_sim_run_task_writes_artifact(tmp_path) -> None:
    cfg = {
        "sim": {
            "name": "dummy",
            "time": {"start": 0.0, "stop": 1.0, "steps": 2},
            "species": ["A", "B", "C"],
        }
    }
    store = ArtifactStore(tmp_path / "artifacts")
    task = get("task", "sim.run")
    backend = resolve_backend("dummy")

    result = task(cfg, store=store)

    assert backend is not None
    assert result.manifest.kind == "runs"
    assert result.path.exists()
    state_path = result.path / "state.zarr" / "dataset.json"
    assert state_path.exists()
    payload = json.loads(state_path.read_text(encoding="utf-8"))
    assert payload["coords"]["species"]["data"] == ["A", "B", "C"]


def test_sim_run_reuses_existing_artifact_before_backend_run(tmp_path) -> None:
    class CountingDummyBackend(DummyBackend):
        name = "dummy_cache_counter"

        def __init__(self) -> None:
            self.calls = 0

        def run(self, cfg):
            self.calls += 1
            return super().run(cfg)

    backend = CountingDummyBackend()
    register("backend", backend.name, backend, overwrite=True)
    store = ArtifactStore(tmp_path / "artifacts")
    task = get("task", "sim.run")
    cfg = {
        "sim": {
            "name": backend.name,
            "time": {"start": 0.0, "stop": 1.0, "steps": 2},
            "species": ["A", "B"],
        }
    }

    first = task(cfg, store=store)
    second = task(cfg, store=store)

    assert first.reused is False
    assert second.reused is True
    assert first.manifest.id == second.manifest.id
    assert backend.calls == 1
