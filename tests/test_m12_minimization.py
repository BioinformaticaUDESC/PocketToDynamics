#!/usr/bin/env python3

import json
from pathlib import Path


JOB = Path(
    "/mnt/labio-fast-disk/"
    "PocketToDynamics/dynamics/"
    "minimized_fixture/"
    "MDTEST_a5b7f2f48479"
)

manifest = (
    JOB /
    "minimization_manifest.json"
)

assert manifest.exists()

data = json.loads(
    manifest.read_text()
)

assert data["status"] == "PASS"

assert data["total_atoms"] == 64587

assert abs(
    data["net_charge_after_ions"]
) < 1e-4

assert (
    data["minimization"][
        "converged"
    ]
    is True
)

assert (
    data["minimization"][
        "potential_energy_kj_mol"
    ]
    < 0
)

for name in [
    "solvated.prmtop",
    "solvated.inpcrd",
    "system_solvated.top",
    "system_solvated.gro",
    "em.tpr",
    "em.gro",
    "em.edr",
    "em.log",
]:

    path = JOB / name

    assert path.exists(), path
    assert path.stat().st_size > 0


text = (
    JOB /
    "em.log"
).read_text(
    errors="replace"
)

assert (
    "converged to Fmax < 1000"
    in text
)

assert (
    "Potential Energy"
    in text
)

assert (
    "Maximum force"
    in text
)


print(
    "PASS_M12_MINIMIZATION_QC"
)

print(
    "Atoms:",
    data["total_atoms"],
)

print(
    "Charge:",
    data[
        "net_charge_after_ions"
    ],
)

print(
    "Potential energy:",
    data["minimization"][
        "potential_energy_kj_mol"
    ],
)

print(
    "Converged:",
    data["minimization"][
        "converged"
    ],
)
