#!/usr/bin/env python3

import json
from pathlib import Path

ROOT = Path.home() / "PocketToDynamics"

JOB = Path(
    "/mnt/labio-fast-disk/"
    "PocketToDynamics/dynamics/"
    "parameterized_fixture/"
    "MDTEST_a5b7f2f48479"
)

manifest = (
    JOB /
    "parameterization_manifest.json"
)

required = [
    "ligand_gaff2.mol2",
    "ligand.frcmod",
    "complex.prmtop",
    "complex.inpcrd",
    "complex_amber.pdb",
    "system.top",
    "system.gro",
    "receptor_clean.pdb",
]

assert manifest.exists(), manifest

data = json.loads(
    manifest.read_text()
)

assert data["status"] == "PASS"

assert (
    data["protein_forcefield"]
    == "ff19SB"
)

assert (
    data["ligand_forcefield"]
    == "gaff2"
)

assert (
    data["charge_method"]
    == "bcc"
)

assert (
    data["receptor_preparation"][
        "strip_crystallographic_waters"
    ]
    is True
)

assert (
    data["receptor_preparation"][
        "removed_water_residues"
    ]
    == 147
)

for name in required:

    path = JOB / name

    assert path.exists(), path
    assert path.stat().st_size > 0


gro = (
    JOB /
    "system.gro"
)

lines = gro.read_text(
    errors="replace"
).splitlines()

assert len(lines) > 3

natoms = int(
    lines[1].strip()
)

assert natoms == 6691


print(
    "PASS_M11_PARAMETERIZATION_QC"
)

print(
    "Atoms:",
    natoms,
)

print(
    "Waters removed:",
    data["receptor_preparation"][
        "removed_water_residues"
    ],
)

print(
    "Protein FF:",
    data[
        "protein_forcefield"
    ],
)

print(
    "Ligand FF:",
    data[
        "ligand_forcefield"
    ],
)
