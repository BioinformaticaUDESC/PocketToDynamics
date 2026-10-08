#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import yaml

from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors


ROOT = Path.home() / "PocketToDynamics"

CONFIG = (
    ROOT /
    "config/chemoinformatics.yaml"
)

RANKING = (
    ROOT /
    "work/ranking/candidate_ranking.tsv"
)

LIGAND_DIR = (
    ROOT /
    "work/normalized/ligands"
)

OUT_DIR = (
    ROOT /
    "work/chemoinformatics"
)

DESCRIPTOR_TABLE = (
    OUT_DIR /
    "ligand_descriptors.tsv"
)

INTEGRATED_TABLE = (
    OUT_DIR /
    "candidate_chemoinformatics.tsv"
)

MODULE_LOG = (
    ROOT /
    "logs/M08_chemoinformatics.log"
)


SUPPORTED = {
    ".sdf",
    ".mol",
    ".mol2",
}


def load_tsv(path: Path):

    with path.open() as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


def load_molecule(path: Path):

    suffix = path.suffix.lower()

    mol = None

    if suffix == ".sdf":

        supplier = Chem.SDMolSupplier(
            str(path),
            removeHs=False,
            sanitize=True,
        )

        for candidate in supplier:

            if candidate is not None:
                mol = candidate
                break

    elif suffix == ".mol":

        mol = Chem.MolFromMolFile(
            str(path),
            removeHs=False,
            sanitize=True,
        )

    elif suffix == ".mol2":

        mol = Chem.MolFromMol2File(
            str(path),
            removeHs=False,
            sanitize=True,
        )

    return mol


def locate_ligand(
    ligand_id: str,
) -> Path | None:

    candidates = []

    for path in LIGAND_DIR.iterdir():

        if not path.is_file():
            continue

        if path.suffix.lower() not in SUPPORTED:
            continue

        if path.stem == ligand_id:
            candidates.append(path)

    if not candidates:
        return None

    candidates.sort()

    return candidates[0]


def molecular_formula(
    mol: Chem.Mol,
) -> str:

    return (
        rdMolDescriptors.CalcMolFormula(
            mol
        )
    )


def formal_charge(
    mol: Chem.Mol,
) -> int:

    return sum(
        atom.GetFormalCharge()
        for atom in mol.GetAtoms()
    )


def descriptors(
    ligand_id: str,
    path: Path,
    mol: Chem.Mol,
):

    heavy_atoms = int(
        mol.GetNumHeavyAtoms()
    )

    return {
        "ligand_id":
            ligand_id,

        "source_file":
            str(
                path.relative_to(ROOT)
            ),

        "molecular_formula":
            molecular_formula(mol),

        "molecular_weight":
            float(
                Descriptors.MolWt(mol)
            ),

        "exact_molecular_weight":
            float(
                Descriptors.ExactMolWt(
                    mol
                )
            ),

        "heavy_atoms":
            heavy_atoms,

        "hbond_donors":
            int(
                Lipinski.NumHDonors(
                    mol
                )
            ),

        "hbond_acceptors":
            int(
                Lipinski.NumHAcceptors(
                    mol
                )
            ),

        "clogp":
            float(
                Crippen.MolLogP(mol)
            ),

        "tpsa":
            float(
                rdMolDescriptors.CalcTPSA(
                    mol
                )
            ),

        "rotatable_bonds":
            int(
                Lipinski.NumRotatableBonds(
                    mol
                )
            ),

        "ring_count":
            int(
                Lipinski.RingCount(
                    mol
                )
            ),

        "aromatic_ring_count":
            int(
                Lipinski.NumAromaticRings(
                    mol
                )
            ),

        "fraction_csp3":
            float(
                rdMolDescriptors.CalcFractionCSP3(
                    mol
                )
            ),

        "formal_charge":
            formal_charge(mol),
    }


def fmt(value):

    if value is None:
        return ""

    if isinstance(value, float):

        if math.isnan(value):
            return ""

        return f"{value:.6f}"

    return value


def main():

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not RANKING.exists():

        print(
            "FAIL_M08_RANKING_NOT_FOUND"
        )

        sys.exit(1)

    cfg = yaml.safe_load(
        CONFIG.read_text()
    )

    reject_zero = bool(
        cfg[
            "chemoinformatics"
        ][
            "qc"
        ][
            "reject_zero_heavy_atoms"
        ]
    )

    ranking_rows = load_tsv(
        RANKING
    )

    ligand_ids = sorted(
        {
            row["ligand_id"]
            for row in ranking_rows
        }
    )

    descriptor_rows = []

    failures = []

    descriptor_lookup = {}

    for ligand_id in ligand_ids:

        path = locate_ligand(
            ligand_id
        )

        if path is None:

            failures.append(
                (
                    ligand_id,
                    "LIGAND_NOT_FOUND",
                )
            )

            print(
                f"[M08] {ligand_id} "
                "FAIL ligand not found"
            )

            continue

        mol = load_molecule(path)

        if mol is None:

            failures.append(
                (
                    ligand_id,
                    "RDKIT_PARSE_FAILED",
                )
            )

            print(
                f"[M08] {ligand_id} "
                "FAIL RDKit parse"
            )

            continue

        row = descriptors(
            ligand_id,
            path,
            mol,
        )

        if (
            reject_zero
            and
            row["heavy_atoms"] == 0
        ):

            failures.append(
                (
                    ligand_id,
                    "ZERO_HEAVY_ATOMS",
                )
            )

            print(
                f"[M08] {ligand_id} "
                "FAIL zero heavy atoms"
            )

            continue

        descriptor_rows.append(
            row
        )

        descriptor_lookup[
            ligand_id
        ] = row

        print(
            f"[M08] {ligand_id} "
            f"PASS "
            f"MW={row['molecular_weight']:.3f} "
            f"HA={row['heavy_atoms']} "
            f"cLogP={row['clogp']:.3f}"
        )

    descriptor_fields = [
        "ligand_id",
        "source_file",
        "molecular_formula",
        "molecular_weight",
        "exact_molecular_weight",
        "heavy_atoms",
        "hbond_donors",
        "hbond_acceptors",
        "clogp",
        "tpsa",
        "rotatable_bonds",
        "ring_count",
        "aromatic_ring_count",
        "fraction_csp3",
        "formal_charge",
    ]

    with DESCRIPTOR_TABLE.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=descriptor_fields,
            delimiter="\t",
        )

        writer.writeheader()

        for row in descriptor_rows:

            writer.writerow(
                {
                    key: fmt(
                        row.get(key)
                    )
                    for key
                    in descriptor_fields
                }
            )

    integrated_rows = []

    for rank_row in ranking_rows:

        ligand_id = (
            rank_row["ligand_id"]
        )

        desc = (
            descriptor_lookup.get(
                ligand_id
            )
        )

        if desc is None:
            continue

        best_score = float(
            rank_row["best_score"]
        )

        heavy_atoms = int(
            desc["heavy_atoms"]
        )

        docking_le = (
            -best_score / heavy_atoms
            if heavy_atoms > 0
            else math.nan
        )

        merged = {
            **rank_row,
            **desc,

            "docking_ligand_efficiency":
                docking_le,
        }

        integrated_rows.append(
            merged
        )

    integrated_fields = [
        "integrated_rank",
        "target_id",
        "ligand_id",

        "n_pockets",
        "n_hits",
        "hit_fraction",

        "best_score",
        "mean_score",
        "median_score",
        "score_sd",
        "worst_score",

        "best_pocket",
        "pocket_selectivity_gap",

        "molecular_formula",
        "molecular_weight",
        "exact_molecular_weight",
        "heavy_atoms",

        "hbond_donors",
        "hbond_acceptors",

        "clogp",
        "tpsa",
        "rotatable_bonds",

        "ring_count",
        "aromatic_ring_count",
        "fraction_csp3",

        "formal_charge",

        "docking_ligand_efficiency",

        "source_file",
    ]

    with INTEGRATED_TABLE.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=integrated_fields,
            delimiter="\t",
        )

        writer.writeheader()

        for row in integrated_rows:

            writer.writerow(
                {
                    key: fmt(
                        row.get(key)
                    )
                    for key
                    in integrated_fields
                }
            )

    if failures:

        status = (
            "FAIL_M08_CHEMOINFORMATICS"
        )

    elif not descriptor_rows:

        status = (
            "FAIL_M08_NO_VALID_LIGANDS"
        )

    else:

        status = (
            "PASS_M08_CHEMOINFORMATICS"
        )

    MODULE_LOG.write_text(
        "\n".join(
            [
                status,
                (
                    "Ligands requested: "
                    f"{len(ligand_ids)}"
                ),
                (
                    "Ligands parsed: "
                    f"{len(descriptor_rows)}"
                ),
                (
                    "Integrated rows: "
                    f"{len(integrated_rows)}"
                ),
                (
                    "Failures: "
                    f"{len(failures)}"
                ),
                (
                    "Descriptor table: "
                    f"{DESCRIPTOR_TABLE}"
                ),
                (
                    "Integrated table: "
                    f"{INTEGRATED_TABLE}"
                ),
            ]
        )
        + "\n"
    )

    print()
    print(status)

    print(
        "Ligands requested:",
        len(ligand_ids),
    )

    print(
        "Ligands parsed:",
        len(descriptor_rows),
    )

    print(
        "Integrated rows:",
        len(integrated_rows),
    )

    if failures:

        print()
        print("Failures:")

        for ligand_id, reason in failures:
            print(
                ligand_id,
                reason,
            )

        sys.exit(1)


if __name__ == "__main__":
    main()
