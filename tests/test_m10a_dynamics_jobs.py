#!/usr/bin/env python3

import csv
import json
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

TABLE = (
    ROOT /
    "work/dynamics/jobs/"
    "dynamics_job_manifest.tsv"
)

SUMMARY = (
    ROOT /
    "work/dynamics/jobs/"
    "dynamics_job_summary.json"
)


assert TABLE.exists()
assert SUMMARY.exists()


with TABLE.open() as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


summary = json.loads(
    SUMMARY.read_text()
)


assert len(rows) == 0

assert (
    summary[
        "promoted_candidates"
    ]
    == 0
)

assert (
    summary[
        "dynamics_jobs"
    ]
    == 0
)

assert (
    summary["status"]
    ==
    "SKIP_M10_NO_PROMOTED_CANDIDATES"
)


print(
    "PASS_M10A_DYNAMICS_JOB_QC"
)

print(
    "Promoted candidates:",
    summary[
        "promoted_candidates"
    ],
)

print(
    "Dynamics jobs:",
    summary[
        "dynamics_jobs"
    ],
)
