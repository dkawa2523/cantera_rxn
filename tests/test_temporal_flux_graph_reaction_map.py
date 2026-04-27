import numpy as np
import pytest


def test_temporal_flux_reaction_map_prefers_index_alignment_when_lengths_match():
    from rxn_platform.graphs import temporal_flux_graph as tfg

    class Stoich:
        def __init__(self):
            # 2 species x 2 reactions
            self.matrix = np.array([[-1.0, 0.0], [1.0, -2.0]], dtype=float)
            self.species = ["A", "B"]
            # Typical stoich IDs (synthetic)
            self.reaction_ids = ["R1", "R2"]
            # Typical Cantera labels (equations)
            self.reaction_equations = ["A => B", "B => 2 A"]

    stoich = Stoich()
    # Run reaction axis uses equations, not synthetic IDs, but ordering is aligned.
    reaction_axis = list(stoich.reaction_equations)
    species_axis = list(stoich.species)

    reaction_map, meta = tfg._build_reaction_map_from_stoich(  # type: ignore[attr-defined]
        stoich, reaction_axis, species_axis
    )
    assert meta["method"] == "index"
    assert len(reaction_map) == 2
    assert all(reaction_map)


def test_temporal_flux_reaction_map_can_fallback_to_equation_matching():
    from rxn_platform.graphs import temporal_flux_graph as tfg

    class Stoich:
        def __init__(self):
            self.matrix = np.array([[-1.0, 0.0], [1.0, -2.0]], dtype=float)
            self.species = ["A", "B"]
            self.reaction_ids = ["R1", "R2"]
            self.reaction_equations = ["eq0", "eq1"]

    stoich = Stoich()
    # Mismatched length: disables index alignment, but equations still match.
    reaction_axis = ["eq0", "eq1", "eq2"]
    species_axis = list(stoich.species)

    reaction_map, meta = tfg._build_reaction_map_from_stoich(  # type: ignore[attr-defined]
        stoich, reaction_axis, species_axis
    )
    assert meta["method"] == "equation"
    assert len(reaction_map) == 3
    # First two reactions should have species entries; the 3rd is absent in stoich.
    assert reaction_map[0]
    assert reaction_map[1]
    assert reaction_map[2] == []


def test_temporal_flux_rejects_missing_stoich_source(tmp_path):
    from rxn_platform.errors import ConfigError
    from rxn_platform.graphs import temporal_flux_graph as tfg

    payload = {
        "coords": {
            "time": {"dims": ["time"], "data": [0.0, 1.0]},
            "species": {"dims": ["species"], "data": ["A", "B"]},
            "reaction": {"dims": ["reaction"], "data": ["R1"]},
        },
        "data_vars": {
            "rop_net": {
                "dims": ["time", "reaction"],
                "data": [[1.0], [2.0]],
            }
        },
        "attrs": {},
    }

    with pytest.raises(ConfigError, match="require mechanism"):
        tfg.build_temporal_flux_graph(
            run_payloads=[payload],
            run_ids=["run-1"],
            mechanism=None,
            phase=None,
            windowing_cfg={"type": "fixed", "count": 1, "base": 10.0},
            aggregation_cfg={"mean": 1.0, "p95": 0.0},
            sparsify_cfg={"quantile": None, "min_value": None},
            rop_cfg={"var": "rop_net", "use_abs": True},
            reaction_map_cfg=None,
            base_dir=tmp_path,
        )

