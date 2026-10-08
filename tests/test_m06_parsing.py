#!/usr/bin/env python3

import csv
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

POSES = (
    ROOT /
    "work/parsed/docking_poses.tsv"
)

TASKS = (
    ROOT /
    "work/parsed/task_summary.tsv"
)


assert POSES.exists()
assert TASKS.exists()


with POSES.open() as f:

    poses = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


with TASKS.open() as f:

    tasks = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )


assert len(tasks) == 1

assert len(poses) == 9

task = tasks[0]

assert task["task_id"] == (
    "TASK_a5b7f2f48479"
)

assert int(
    task["n_modes"]
) == 9


scores = [
    float(row["affinity_kcal_mol"])
    for row in poses
]

best = min(scores)

assert abs(
    float(task["best_score"])
    - best
) < 1e-6


best_modes = [
    row
    for row in poses
    if row["is_best_mode"] == "1"
]

assert len(best_modes) >= 1


print("PASS_M06_PARSING_QC")
print("Tasks:", len(tasks))
print("Poses:", len(poses))
print("Best score:", best)
