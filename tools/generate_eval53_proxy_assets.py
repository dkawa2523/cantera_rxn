"""Generate transparent Cantera-readable proxy assets for eval53 large reruns.

The original eval53 large Cantera mechanisms are not present in this checkout.
This script creates local proxy assets that let the reduction-sweep workflow run
end-to-end without pretending to recover the original benchmark chemistry.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import cantera as ct


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "benchmarks" / "assets" / "eval53_large"


def _species(name: str, composition: dict[str, float]) -> ct.Species:
    species = ct.Species(name, composition)
    species.thermo = ct.ConstantCp(
        200.0,
        5000.0,
        ct.one_atm,
        [298.15, 0.0, 0.0, 33_000.0],
    )
    return species


def _reaction(equation: str, index: int) -> ct.Reaction:
    rate = ct.ArrheniusRate(2.5e4 * (1.0 + 0.02 * index), 0.0, 12_000.0 + 250.0 * index)
    return ct.Reaction(equation=equation, rate=rate)


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _generate_ac() -> None:
    source = ROOT / "benchmarks" / "assets" / "mechanisms" / "gri30.yaml"
    target = ASSETS / "mechanisms" / "ac_large.yaml"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
    _write_csv(
        ASSETS / "conditions" / "ac.csv",
        [
            {
                "case_id": "dilute_N2",
                "T": "1050.0",
                "P": "1333.2236842105265",
                "t_end": "0.05",
                "gas_X": "C2H2:0.02,O2:0.05,H2:0.03,N2:0.90",
            }
        ],
    )


def _generate_sif4() -> None:
    specs = [
        _species("H", {"H": 1}),
        _species("H2", {"H": 2}),
        _species("O", {"O": 1}),
        _species("O2", {"O": 2}),
        _species("OH", {"O": 1, "H": 1}),
        _species("HO2", {"H": 1, "O": 2}),
        _species("H2O", {"H": 2, "O": 1}),
        _species("H2O2", {"H": 2, "O": 2}),
        _species("N", {"N": 1}),
        _species("N2", {"N": 2}),
        _species("NH", {"N": 1, "H": 1}),
        _species("NH2", {"N": 1, "H": 2}),
        _species("NH3", {"N": 1, "H": 3}),
        _species("NO", {"N": 1, "O": 1}),
        _species("HNO", {"H": 1, "N": 1, "O": 1}),
        _species("F", {"F": 1}),
        _species("F2", {"F": 2}),
        _species("HF", {"H": 1, "F": 1}),
        _species("Si", {"Si": 1}),
        _species("SiF", {"Si": 1, "F": 1}),
        _species("SiF2", {"Si": 1, "F": 2}),
        _species("SiF3", {"Si": 1, "F": 3}),
        _species("SiF4", {"Si": 1, "F": 4}),
        _species("SiO", {"Si": 1, "O": 1}),
        _species("SiO2", {"Si": 1, "O": 2}),
        _species("SiFO", {"Si": 1, "F": 1, "O": 1}),
    ]
    equations = [
        "SiF4 + H <=> SiF3 + HF",
        "SiF3 + H <=> SiF2 + HF",
        "SiF2 + H <=> SiF + HF",
        "SiF + H <=> Si + HF",
        "NH3 + H <=> NH2 + H2",
        "NH2 + H <=> NH + H2",
        "NH + H <=> N + H2",
        "O2 + H <=> O + OH",
        "OH + H2 <=> H2O + H",
        "O + H2 <=> OH + H",
        "NH2 + O <=> HNO + H",
        "HNO + H <=> NO + H2",
        "NO + N <=> N2 + O",
        "Si + O2 <=> SiO + O",
        "SiO + O <=> SiO2",
        "SiF2 + O <=> SiFO + F",
        "F + H2 <=> HF + H",
        "F2 + H <=> HF + F",
        "SiF3 + F <=> SiF4",
        "SiF2 + F <=> SiF3",
        "SiF + F <=> SiF2",
        "Si + F <=> SiF",
        "NH + O <=> NO + H",
        "N + O2 <=> NO + O",
        "HO2 + H <=> 2 OH",
        "H + O2 <=> HO2",
        "H2O2 + H <=> H2O + OH",
        "2 OH <=> H2O2",
        "SiF4 + NH3 <=> SiF3 + NH2 + HF",
        "SiF3 + NH2 <=> SiF2 + NH + HF",
        "SiF2 + NH <=> SiF + N + HF",
        "SiF + O <=> SiO + F",
        "SiFO + H <=> SiO + HF",
        "NH3 + OH <=> NH2 + H2O",
        "NH2 + OH <=> NH + H2O",
        "N + N <=> N2",
        "H + H <=> H2",
        "O + O <=> O2",
        "F + F <=> F2",
        "SiO + OH <=> SiO2 + H",
        "SiF2 + OH <=> SiFO + HF",
    ]
    gas = ct.Solution(
        thermo="ideal-gas",
        kinetics="gas",
        species=specs,
        reactions=[_reaction(equation, index) for index, equation in enumerate(equations)],
        name="gas",
    )
    gas.TPX = 1713.0, 266.6447368421053, "SiF4:0.08,NH3:0.16,O2:0.02,N2:0.74"
    target = ASSETS / "mechanisms" / "sif4_large.yaml"
    target.parent.mkdir(parents=True, exist_ok=True)
    gas.write_yaml(str(target), header=True)
    _write_csv(
        ASSETS / "conditions" / "sif4.csv",
        [
            {
                "case_id": "add_O2",
                "T": "1713.0",
                "P": "266.6447368421053",
                "t_end": "0.002",
                "gas_X": "SiF4:0.08,NH3:0.16,O2:0.02,N2:0.74",
            }
        ],
    )


def _write_readme() -> None:
    payload = {
        "kind": "generated_proxy_assets",
        "note": (
            "These assets are Cantera-readable proxies generated because the "
            "original eval53 large mechanisms are not present in this checkout."
        ),
        "ac_large.yaml": "copy of local benchmarks/assets/mechanisms/gri30.yaml",
        "sif4_large.yaml": "synthetic constant-cp SiF4/NH3/O2 gas-phase proxy",
        "conditions": ["ac:dilute_N2", "sif4:add_O2"],
    }
    (ASSETS / "proxy_assets.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (ASSETS / "README.md").write_text(
        "# Eval53 Large Proxy Assets\n\n"
        "The original eval53 large Cantera mechanisms are not present in this checkout.\n"
        "These files were generated as transparent Cantera-readable proxies so the\n"
        "learnCK compression-sweep workflow can run end-to-end.\n\n"
        "- `mechanisms/ac_large.yaml`: local GRI30 gas mechanism copy used as an acetylene proxy.\n"
        "- `mechanisms/sif4_large.yaml`: synthetic constant-cp SiF4/NH3/O2 gas-phase proxy.\n"
        "- `conditions/ac.csv`: contains `dilute_N2`.\n"
        "- `conditions/sif4.csv`: contains `add_O2`.\n\n"
        "Do not treat these files as recovered original eval53 large mechanisms.\n",
        encoding="utf-8",
    )


def main() -> None:
    _generate_ac()
    _generate_sif4()
    _write_readme()
    print(f"wrote proxy assets under {ASSETS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
