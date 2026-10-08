#!/usr/bin/env python3

import json
from pathlib import Path


JOB = Path(
    "/mnt/labio-fast-disk/"
    "PocketToDynamics/dynamics/"
    "equilibrated_fixture/"
    "MDTEST_a5b7f2f48479"
)

manifest = (
    JOB /
    "equilibration_manifest.json"
)

assert manifest.exists(), manifest

data = json.loads(
    manifest.read_text()
)

assert data["status"] == "PASS"

target_t = float(
    data["temperature_target_k"]
)

nvt_t = float(
    data["nvt"]["temperature"]["mean"]
)

npt_t = float(
    data["npt"]["temperature"]["mean"]
)

density = float(
    data["npt"]["density"]["mean"]
)

assert abs(
    nvt_t - target_t
) < 10.0

assert abs(
    npt_t - target_t
) < 10.0

assert (
    900.0
    <= density
    <= 1100.0
)


required = [
    "nvt.gro",
    "nvt.cpt",
    "nvt.edr",
    "nvt.log",
    "npt.gro",
    "npt.cpt",
    "npt.edr",
    "npt.log",
    "system.top",
]

for name in required:

    path = JOB / name

    assert path.exists(), path
    assert path.stat().st_size > 0


warning_patterns = [
    "LINCS WARNING",
    "constraint warning",
    "nan occurred",
    "fatal error",
    "segmentation fault",
    "blowing up",
]


import re

constraint_deviations = []

for filename in [
    "nvt.log",
    "npt.log",
]:

    log_text = (
        JOB /
        filename
    ).read_text(
        errors="replace"
    )

    matches = re.findall(
        r"RMS relative constraint deviation "
        r"after constraining:\s*"
        r"([-+0-9.eE]+)",
        log_text,
    )

    for value in matches:
        constraint_deviations.append(
            float(value)
        )

assert constraint_deviations, (
    "No constraint-deviation diagnostics found"
)

assert max(constraint_deviations) < 1e-3, (
    "Constraint deviation too large: "
    f"{max(constraint_deviations)}"
)

for filename in [
    "nvt.log",
    "npt.log",
    "mdrun_nvt.stderr",
    "mdrun_npt.stderr",
]:

    text = (
        JOB /
        filename
    ).read_text(
        errors="replace"
    ).lower()

    for pattern in warning_patterns:

        assert (
            pattern.lower()
            not in text
        ), (
            f"{pattern} found "
            f"in {filename}"
        )


print(
    "PASS_M13_EQUILIBRATION_QC"
)

print(
    "Maximum RMS constraint deviation:",
    max(constraint_deviations),
)

print(
    "NVT temperature:",
    nvt_t,
)

print(
    "NPT temperature:",
    npt_t,
)

print(
    "NPT density:",
    density,
)

print(
    "NPT pressure:",
    data["npt"]["pressure"]["mean"],
)
