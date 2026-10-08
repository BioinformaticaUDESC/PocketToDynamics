#!/usr/bin/env python3

import csv
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

TABLE = (
    ROOT /
    "work/ranking/candidate_ranking.tsv"
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


assert (
    row["target_id"]
    == "test_receptor"
)

assert (
    row["ligand_id"]
    == "test_ligand"
)

assert int(
    row["n_pockets"]
) == 1

assert int(
    row["n_hits"]
) == 0

assert abs(
    float(
        row["best_score"]
    )
    - (-2.399)
) < 1e-6

assert (
    row["best_pocket"]
    == "P0001"
)

assert int(
    row["integrated_rank"]
) == 1


print(
    "PASS_M07_RANKING_QC"
)

print(
    "Candidates:",
    len(rows),
)

print(
    "Best score:",
    row["best_score"],
)

print(
    "Best pocket:",
    row["best_pocket"],
)
