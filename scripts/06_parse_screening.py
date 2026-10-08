#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

RAW_ROOT = ROOT / "work" / "screening" / "raw"
PARSED_ROOT = ROOT / "work" / "parsed"

POSES_TSV = PARSED_ROOT / "docking_poses.tsv"
TASKS_TSV = PARSED_ROOT / "task_summary.tsv"

MODULE_LOG = ROOT / "logs" / "M06_parsing.log"


VINA_ROW = re.compile(
    r"^\s*"
    r"(?P<mode>\d+)"
    r"\s+"
    r"(?P<affinity>-?\d+(?:\.\d+)?)"
    r"\s+"
    r"(?P<rmsd_lb>-?\d+(?:\.\d+)?)"
    r"\s+"
    r"(?P<rmsd_ub>-?\d+(?:\.\d+)?)"
    r"\s*$"
)


def mean(values):
    return (
        statistics.mean(values)
        if values
        else math.nan
    )


def median(values):
    return (
        statistics.median(values)
        if values
        else math.nan
    )


def stdev(values):
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def parse_vina_log(path: Path):

    rows = []

    for line in path.read_text(
        errors="replace"
    ).splitlines():

        match = VINA_ROW.match(line)

        if not match:
            continue

        rows.append(
            {
                "mode": int(
                    match.group("mode")
                ),
                "affinity_kcal_mol": float(
                    match.group("affinity")
                ),
                "rmsd_lb": float(
                    match.group("rmsd_lb")
                ),
                "rmsd_ub": float(
                    match.group("rmsd_ub")
                ),
            }
        )

    return rows


def fmt(value):

    if isinstance(value, float):
        return f"{value:.6f}"

    return value


def main():

    PARSED_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    manifests = sorted(
        RAW_ROOT.glob(
            "*/task_manifest.json"
        )
    )

    pose_rows = []
    task_rows = []

    pass_tasks = 0
    fail_tasks = 0

    for manifest_path in manifests:

        manifest = json.loads(
            manifest_path.read_text()
        )

        task_id = manifest["task_id"]

        if manifest.get("status") != "PASS":

            print(
                f"[M06] {task_id} "
                f"SKIP_NONPASS_TASK"
            )

            continue

        task_dir = manifest_path.parent

        vina_log = (
            task_dir /
            "vina.log"
        )

        if not vina_log.exists():

            print(
                f"[M06] {task_id} "
                f"FAIL_NO_LOG"
            )

            fail_tasks += 1
            continue

        parsed = parse_vina_log(
            vina_log
        )

        if not parsed:

            print(
                f"[M06] {task_id} "
                f"FAIL_NO_POSES_PARSED"
            )

            fail_tasks += 1
            continue

        affinities = [
            row["affinity_kcal_mol"]
            for row in parsed
        ]

        best = min(affinities)

        for row in parsed:

            pose_rows.append(
                {
                    "task_id": task_id,
                    "target_id": manifest[
                        "target_id"
                    ],
                    "pocket_id": manifest[
                        "pocket_id"
                    ],
                    "ligand_id": manifest[
                        "ligand_id"
                    ],
                    "engine": manifest[
                        "engine"
                    ],
                    "engine_version": manifest[
                        "engine_version"
                    ],
                    "mode": row["mode"],
                    "affinity_kcal_mol":
                        row[
                            "affinity_kcal_mol"
                        ],
                    "rmsd_lb":
                        row["rmsd_lb"],
                    "rmsd_ub":
                        row["rmsd_ub"],
                    "is_best_mode":
                        int(
                            row[
                                "affinity_kcal_mol"
                            ] == best
                        ),
                }
            )

        task_rows.append(
            {
                "task_id": task_id,
                "target_id":
                    manifest["target_id"],
                "pocket_id":
                    manifest["pocket_id"],
                "ligand_id":
                    manifest["ligand_id"],
                "engine":
                    manifest["engine"],
                "n_modes":
                    len(parsed),
                "best_score":
                    min(affinities),
                "mean_score":
                    mean(affinities),
                "median_score":
                    median(affinities),
                "score_sd":
                    stdev(affinities),
                "worst_score":
                    max(affinities),
                "runtime_seconds":
                    manifest.get(
                        "runtime_seconds",
                        "",
                    ),
            }
        )

        pass_tasks += 1

        print(
            f"[M06] {task_id} "
            f"PASS "
            f"modes={len(parsed)} "
            f"best={best:.3f}"
        )

    pose_fields = [
        "task_id",
        "target_id",
        "pocket_id",
        "ligand_id",
        "engine",
        "engine_version",
        "mode",
        "affinity_kcal_mol",
        "rmsd_lb",
        "rmsd_ub",
        "is_best_mode",
    ]

    task_fields = [
        "task_id",
        "target_id",
        "pocket_id",
        "ligand_id",
        "engine",
        "n_modes",
        "best_score",
        "mean_score",
        "median_score",
        "score_sd",
        "worst_score",
        "runtime_seconds",
    ]

    with POSES_TSV.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=pose_fields,
            delimiter="\t",
        )

        writer.writeheader()

        for row in pose_rows:
            writer.writerow(
                {
                    k: fmt(v)
                    for k, v in row.items()
                }
            )

    with TASKS_TSV.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=task_fields,
            delimiter="\t",
        )

        writer.writeheader()

        for row in task_rows:
            writer.writerow(
                {
                    k: fmt(v)
                    for k, v in row.items()
                }
            )

    if not manifests:

        status = (
            "SKIP_M06_NO_SCREENING_RESULTS"
        )

    elif fail_tasks > 0:

        status = (
            "FAIL_M06_PARSING"
        )

    elif pass_tasks == 0:

        status = (
            "FAIL_M06_NO_VALID_TASKS"
        )

    else:

        status = (
            "PASS_M06_PARSING"
        )

    MODULE_LOG.write_text(
        "\n".join(
            [
                status,
                (
                    "Task manifests: "
                    f"{len(manifests)}"
                ),
                (
                    "Parsed tasks: "
                    f"{pass_tasks}"
                ),
                (
                    "Failed tasks: "
                    f"{fail_tasks}"
                ),
                (
                    "Parsed poses: "
                    f"{len(pose_rows)}"
                ),
                (
                    f"Poses table: "
                    f"{POSES_TSV}"
                ),
                (
                    f"Task table: "
                    f"{TASKS_TSV}"
                ),
            ]
        )
        + "\n"
    )

    print()
    print(status)
    print(
        "Task manifests:",
        len(manifests),
    )
    print(
        "Parsed tasks:",
        pass_tasks,
    )
    print(
        "Parsed poses:",
        len(pose_rows),
    )

    if status.startswith("FAIL_"):
        sys.exit(1)


if __name__ == "__main__":
    main()
