#!/usr/bin/env python3

import csv
import json
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

BASE = (
    ROOT /
    "work/dynamics/test_fixture"
)

TABLE = (
    BASE /
    "pose_extraction_manifest.tsv"
)

JOB = (
    BASE /
    "MDTEST_a5b7f2f48479"
)


assert TABLE.exists()

with TABLE.open() as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


assert len(rows) == 1

row = rows[0]

assert row["status"] == "PASS"

assert row["selected_mode"] == "1"

assert (
    row["screening_task_id"]
    == "TASK_a5b7f2f48479"
)


required = [
    "selected_pose.pdbqt",
    "selected_pose.mol2",
    "selected_pose.pdb",
    "receptor.pdb",
    "pose_manifest.json",
]


for name in required:

    path = JOB / name

    assert path.exists(), path
    assert path.stat().st_size > 0


text = (
    JOB /
    "selected_pose.pdbqt"
).read_text(
    errors="replace"
)


assert "MODEL 1" in text
assert "MODEL 2" not in text


manifest = json.loads(
    (
        JOB /
        "pose_manifest.json"
    ).read_text()
)


assert manifest["status"] == "PASS"

assert (
    manifest["selected_mode"]
    == 1
)


print(
    "PASS_M10C_POSE_EXTRACTION_QC"
)

print(
    "Job:",
    row["dynamics_job_id"],
)

print(
    "Mode:",
    row["selected_mode"],
)

print(
    "Score:",
    row["selected_score"],
)
