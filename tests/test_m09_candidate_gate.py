#!/usr/bin/env python3

import csv
import json
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

ALL = (
    ROOT /
    "work/candidates/"
    "candidate_gate.tsv"
)

PROMOTED = (
    ROOT /
    "work/candidates/"
    "promoted_candidates.tsv"
)

SUMMARY = (
    ROOT /
    "work/candidates/"
    "candidate_gate_summary.json"
)


assert ALL.exists()
assert PROMOTED.exists()
assert SUMMARY.exists()


with ALL.open() as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


with PROMOTED.open() as f:

    promoted = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


summary = json.loads(
    SUMMARY.read_text()
)


assert len(rows) == 1

row = rows[0]


assert (
    row["target_id"]
    == "test_receptor"
)

assert (
    row["ligand_id"]
    == "test_ligand"
)

assert row["gate_pass"] == "0"

assert len(promoted) == 0

assert (
    summary["evaluated"]
    == 1
)

assert (
    summary["promoted"]
    == 0
)

assert (
    summary["rejected"]
    == 1
)


print(
    "PASS_M09_CANDIDATE_GATE_QC"
)

print(
    "Evaluated:",
    summary["evaluated"],
)

print(
    "Promoted:",
    summary["promoted"],
)

print(
    "Rejected:",
    summary["rejected"],
)

print(
    "Failed criteria:",
    row["failed_criteria"],
)
