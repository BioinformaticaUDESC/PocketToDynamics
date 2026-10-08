#!/usr/bin/env python3

import csv
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

DESC = (
    ROOT /
    "work/chemoinformatics/"
    "ligand_descriptors.tsv"
)

INTEGRATED = (
    ROOT /
    "work/chemoinformatics/"
    "candidate_chemoinformatics.tsv"
)


assert DESC.exists()
assert INTEGRATED.exists()


with DESC.open() as f:

    desc = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


with INTEGRATED.open() as f:

    integrated = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


assert len(desc) == 1
assert len(integrated) == 1


d = desc[0]
r = integrated[0]


assert (
    d["ligand_id"]
    == "test_ligand"
)

assert (
    r["ligand_id"]
    == "test_ligand"
)

assert int(
    d["heavy_atoms"]
) > 0


mw = float(
    d["molecular_weight"]
)

assert mw > 0


le = float(
    r[
        "docking_ligand_efficiency"
    ]
)

expected = (
    -float(
        r["best_score"]
    )
    /
    int(
        r["heavy_atoms"]
    )
)

assert abs(
    le - expected
) < 1e-6


print(
    "PASS_M08_CHEMOINFORMATICS_QC"
)

print(
    "Ligands:",
    len(desc),
)

print(
    "MW:",
    d["molecular_weight"],
)

print(
    "Heavy atoms:",
    d["heavy_atoms"],
)

print(
    "Docking LE:",
    r[
        "docking_ligand_efficiency"
    ],
)
