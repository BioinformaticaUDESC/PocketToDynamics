#!/usr/bin/env python3

from pathlib import Path
import csv
import sys

ROOT = Path.home() / "PocketToDynamics"

POCKETS = (
    ROOT /
    "work/pockets/test_receptor/pockets.tsv"
)

JOBS = (
    ROOT /
    "work/jobs/job_matrix.tsv"
)


def read_tsv(path):
    with path.open(
        newline=""
    ) as handle:
        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def main():

    failures = []

    if not POCKETS.exists():
        failures.append(
            "Missing pockets.tsv"
        )

    if not JOBS.exists():
        failures.append(
            "Missing job_matrix.tsv"
        )

    if failures:
        for failure in failures:
            print(
                "FAIL:",
                failure,
            )

        sys.exit(1)

    pockets = read_tsv(
        POCKETS
    )

    jobs = read_tsv(
        JOBS
    )

    targets = {
        row["target_id"]
        for row in jobs
    }

    ligands = {
        row["ligand_id"]
        for row in jobs
    }

    expected = (
        len(targets)
        * len(pockets)
        * len(ligands)
    )

    if len(jobs) != expected:
        failures.append(
            "Job count mismatch: "
            f"expected={expected}, "
            f"observed={len(jobs)}"
        )

    task_ids = [
        row["task_id"]
        for row in jobs
    ]

    if len(task_ids) != len(
        set(task_ids)
    ):
        failures.append(
            "Duplicate task_id detected"
        )

    if any(
        row["status"] != "PENDING"
        for row in jobs
    ):
        failures.append(
            "Initial jobs must be PENDING"
        )

    for row in jobs:

        for field in (
            "center_x",
            "center_y",
            "center_z",
        ):

            try:
                float(
                    row[field]
                )
            except ValueError:

                failures.append(
                    f"Invalid {field} "
                    f"for {row['task_id']}"
                )

    if failures:

        print(
            "FAIL_INTEGRATION_M01_M04"
        )

        for failure in failures:
            print(
                " -",
                failure,
            )

        sys.exit(1)

    print(
        "PASS_INTEGRATION_M01_M04"
    )

    print(
        f"Targets: {len(targets)}"
    )

    print(
        f"Pockets: {len(pockets)}"
    )

    print(
        f"Ligands: {len(ligands)}"
    )

    print(
        f"Tasks: {len(jobs)}"
    )


if __name__ == "__main__":
    main()
