#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import yaml


ROOT = Path.home() / "PocketToDynamics"

CONFIG = ROOT / "config/dynamics.yaml"

PROMOTED = (
    ROOT /
    "work/candidates/"
    "promoted_candidates.tsv"
)

POSES = (
    ROOT /
    "work/parsed/"
    "docking_poses.tsv"
)

JOB_MATRIX = (
    ROOT /
    "work/jobs/"
    "job_matrix.tsv"
)

SCREENING_ROOT = (
    ROOT /
    "work/screening/raw"
)

OUT_DIR = (
    ROOT /
    "work/dynamics/jobs"
)

OUT_TABLE = (
    OUT_DIR /
    "dynamics_job_manifest.tsv"
)

SUMMARY = (
    OUT_DIR /
    "dynamics_job_summary.json"
)

MODULE_LOG = (
    ROOT /
    "logs/M10A_dynamics_jobs.log"
)


def sha256_text(text: str) -> str:

    return hashlib.sha256(
        text.encode()
    ).hexdigest()


def load_tsv(path: Path):

    if not path.exists():
        return []

    with path.open() as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


def main():

    cfg = yaml.safe_load(
        CONFIG.read_text()
    )

    promoted = load_tsv(
        PROMOTED
    )

    poses = load_tsv(
        POSES
    )

    jobs = load_tsv(
        JOB_MATRIX
    )

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not promoted:

        fields = [
            "dynamics_job_id",
            "target_id",
            "ligand_id",
            "pocket_id",
            "screening_task_id",
            "selected_mode",
            "selected_score",
            "source_pose_file",
            "status",
        ]

        with OUT_TABLE.open(
            "w",
            newline="",
        ) as f:

            writer = csv.DictWriter(
                f,
                fieldnames=fields,
                delimiter="\t",
            )

            writer.writeheader()

        summary = {
            "promoted_candidates": 0,
            "dynamics_jobs": 0,
            "status":
                "SKIP_M10_NO_PROMOTED_CANDIDATES",
        }

        SUMMARY.write_text(
            json.dumps(
                summary,
                indent=2,
            )
            + "\n"
        )

        MODULE_LOG.write_text(
            "\n".join(
                [
                    "SKIP_M10_NO_PROMOTED_CANDIDATES",
                    "Promoted candidates: 0",
                    "Dynamics jobs: 0",
                ]
            )
            + "\n"
        )

        print(
            "SKIP_M10_NO_PROMOTED_CANDIDATES"
        )
        print(
            "Promoted candidates: 0"
        )
        print(
            "Dynamics jobs: 0"
        )

        return

    task_lookup = {}

    for row in jobs:

        key = (
            row["target_id"],
            row["pocket_id"],
            row["ligand_id"],
        )

        task_lookup[key] = row

    pose_lookup = {}

    for row in poses:

        key = row["task_id"]

        if key not in pose_lookup:
            pose_lookup[key] = []

        pose_lookup[key].append(row)

    output_rows = []

    failures = []

    for row in promoted:

        target_id = (
            row["target_id"]
        )

        ligand_id = (
            row["ligand_id"]
        )

        pocket_id = (
            row["best_pocket"]
        )

        key = (
            target_id,
            pocket_id,
            ligand_id,
        )

        job = task_lookup.get(
            key
        )

        if not job:

            failures.append(
                (
                    target_id,
                    ligand_id,
                    "SCREENING_TASK_NOT_FOUND",
                )
            )
            continue

        task_id = (
            job["task_id"]
        )

        task_poses = (
            pose_lookup.get(
                task_id,
                []
            )
        )

        if not task_poses:

            failures.append(
                (
                    target_id,
                    ligand_id,
                    "POSES_NOT_FOUND",
                )
            )
            continue

        best_pose = min(
            task_poses,
            key=lambda x: float(
                x["affinity_kcal_mol"]
            ),
        )

        source_pose_file = (
            SCREENING_ROOT /
            task_id /
            "poses.pdbqt"
        )

        if not source_pose_file.exists():

            failures.append(
                (
                    target_id,
                    ligand_id,
                    "SOURCE_POSE_FILE_NOT_FOUND",
                )
            )
            continue

        identity = (
            f"{target_id}|"
            f"{ligand_id}|"
            f"{pocket_id}|"
            f"{task_id}"
        )

        dynamics_job_id = (
            "MDJOB_"
            + sha256_text(
                identity
            )[:12]
        )

        output_rows.append(
            {
                "dynamics_job_id":
                    dynamics_job_id,

                "target_id":
                    target_id,

                "ligand_id":
                    ligand_id,

                "pocket_id":
                    pocket_id,

                "screening_task_id":
                    task_id,

                "selected_mode":
                    best_pose["mode"],

                "selected_score":
                    best_pose[
                        "affinity_kcal_mol"
                    ],

                "source_pose_file":
                    str(
                        source_pose_file.relative_to(
                            ROOT
                        )
                    ),

                "status":
                    "PENDING",
            }
        )

    fields = [
        "dynamics_job_id",
        "target_id",
        "ligand_id",
        "pocket_id",
        "screening_task_id",
        "selected_mode",
        "selected_score",
        "source_pose_file",
        "status",
    ]

    with OUT_TABLE.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            output_rows
        )

    if failures:

        status = (
            "FAIL_M10_DYNAMICS_JOB_BUILD"
        )

    else:

        status = (
            "PASS_M10_DYNAMICS_JOB_BUILD"
        )

    summary = {
        "promoted_candidates":
            len(promoted),

        "dynamics_jobs":
            len(output_rows),

        "failures":
            len(failures),

        "status":
            status,
    }

    SUMMARY.write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n"
    )

    MODULE_LOG.write_text(
        "\n".join(
            [
                status,
                (
                    "Promoted candidates: "
                    f"{len(promoted)}"
                ),
                (
                    "Dynamics jobs: "
                    f"{len(output_rows)}"
                ),
                (
                    "Failures: "
                    f"{len(failures)}"
                ),
            ]
        )
        + "\n"
    )

    print(status)

    print(
        "Promoted candidates:",
        len(promoted),
    )

    print(
        "Dynamics jobs:",
        len(output_rows),
    )

    print(
        "Failures:",
        len(failures),
    )

    if failures:

        print()
        print("Failure details:")

        for item in failures:
            print(*item)

        sys.exit(1)


if __name__ == "__main__":
    main()
