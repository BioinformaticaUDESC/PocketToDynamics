#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import csv
import json
import hashlib
import sys

ROOT = Path.home() / "PocketToDynamics"

NORMALIZED = ROOT / "work/normalized"
POCKET_ROOT = ROOT / "work/pockets"
JOB_ROOT = ROOT / "work/jobs"
LOG_DIR = ROOT / "logs"

RECEPTOR_MAP = NORMALIZED / "receptor_id_map.tsv"
LIGAND_MAP = NORMALIZED / "ligand_id_map.tsv"

JOB_TABLE = JOB_ROOT / "job_matrix.tsv"
JOB_SUMMARY = JOB_ROOT / "job_summary.json"


def read_tsv(path):
    if not path.exists():
        return []

    with path.open(newline="") as handle:
        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def stable_task_id(target_id, pocket_id, ligand_id):
    raw = f"{target_id}|{pocket_id}|{ligand_id}"

    digest = hashlib.sha256(
        raw.encode("utf-8")
    ).hexdigest()[:12]

    return f"TASK_{digest}"


def load_pockets(target_id):
    path = POCKET_ROOT / target_id / "pockets.tsv"

    if not path.exists():
        return []

    return read_tsv(path)


def main():
    JOB_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    receptors = read_tsv(RECEPTOR_MAP)
    ligands = read_tsv(LIGAND_MAP)

    if not receptors and not ligands:
        status = "SKIP_M04_NO_INPUTS"

        LOG_DIR.joinpath(
            "M04_job_builder.log"
        ).write_text(
            status + "\n"
            "Receptors: 0\n"
            "Ligands: 0\n"
            "Tasks: 0\n"
        )

        print(status)
        print("No normalized receptors or ligands available.")
        return

    if not receptors:
        status = "SKIP_M04_NO_RECEPTORS"

        LOG_DIR.joinpath(
            "M04_job_builder.log"
        ).write_text(
            status + "\n"
            f"Receptors: 0\n"
            f"Ligands: {len(ligands)}\n"
            "Tasks: 0\n"
        )

        print(status)
        return

    if not ligands:
        status = "SKIP_M04_NO_LIGANDS"

        LOG_DIR.joinpath(
            "M04_job_builder.log"
        ).write_text(
            status + "\n"
            f"Receptors: {len(receptors)}\n"
            "Ligands: 0\n"
            "Tasks: 0\n"
        )

        print(status)
        return

    jobs = []
    targets_without_pockets = []

    for receptor in receptors:
        target_id = receptor["internal_id"]

        pockets = load_pockets(target_id)

        if not pockets:
            targets_without_pockets.append(
                target_id
            )
            continue

        for pocket in pockets:
            pocket_id = pocket["pocket_id"]

            for ligand in ligands:
                ligand_id = ligand["internal_id"]

                task_id = stable_task_id(
                    target_id,
                    pocket_id,
                    ligand_id,
                )

                jobs.append({
                    "task_id": task_id,
                    "target_id": target_id,
                    "pocket_id": pocket_id,
                    "ligand_id": ligand_id,

                    "receptor_path":
                        receptor[
                            "normalized_path"
                        ],

                    "ligand_path":
                        ligand[
                            "normalized_path"
                        ],

                    "center_x":
                        pocket.get(
                            "center_x",
                            "",
                        ),

                    "center_y":
                        pocket.get(
                            "center_y",
                            "",
                        ),

                    "center_z":
                        pocket.get(
                            "center_z",
                            "",
                        ),

                    "status": "PENDING",
                    "attempts": "0",
                    "output_path": "",
                    "log_path": "",
                })

    fields = [
        "task_id",
        "target_id",
        "pocket_id",
        "ligand_id",
        "receptor_path",
        "ligand_path",
        "center_x",
        "center_y",
        "center_z",
        "status",
        "attempts",
        "output_path",
        "log_path",
    ]

    with JOB_TABLE.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(jobs)

    unique_targets = sorted({
        j["target_id"]
        for j in jobs
    })

    unique_pockets = sorted({
        (
            j["target_id"],
            j["pocket_id"],
        )
        for j in jobs
    })

    unique_ligands = sorted({
        j["ligand_id"]
        for j in jobs
    })

    if not jobs:
        status = "SKIP_M04_NO_VALID_TASKS"
    elif targets_without_pockets:
        status = "WARN_M04_PARTIAL_JOB_MATRIX"
    else:
        status = "PASS_M04_JOB_MATRIX"

    summary = {
        "module": "M04_JOB_MATRIX_BUILDER",
        "timestamp": datetime.now().isoformat(),
        "status": status,

        "n_receptors_registered":
            len(receptors),

        "n_ligands_registered":
            len(ligands),

        "n_targets_in_jobs":
            len(unique_targets),

        "n_target_pockets":
            len(unique_pockets),

        "n_ligands_in_jobs":
            len(unique_ligands),

        "n_tasks":
            len(jobs),

        "targets_without_pockets":
            targets_without_pockets,

        "job_table":
            str(
                JOB_TABLE.relative_to(ROOT)
            ),
    }

    with JOB_SUMMARY.open("w") as handle:
        json.dump(
            summary,
            handle,
            indent=2,
        )

    log_path = (
        LOG_DIR /
        "M04_job_builder.log"
    )

    with log_path.open("w") as handle:
        handle.write(status + "\n")

        handle.write(
            f"Registered receptors: "
            f"{len(receptors)}\n"
        )

        handle.write(
            f"Registered ligands: "
            f"{len(ligands)}\n"
        )

        handle.write(
            f"Targets represented: "
            f"{len(unique_targets)}\n"
        )

        handle.write(
            f"Target-pockets: "
            f"{len(unique_pockets)}\n"
        )

        handle.write(
            f"Tasks: {len(jobs)}\n"
        )

        handle.write(
            "Targets without pockets: "
            + (
                ",".join(
                    targets_without_pockets
                )
                if targets_without_pockets
                else "none"
            )
            + "\n"
        )

    print(status)
    print(
        f"Registered receptors: "
        f"{len(receptors)}"
    )

    print(
        f"Registered ligands: "
        f"{len(ligands)}"
    )

    print(
        f"Targets represented: "
        f"{len(unique_targets)}"
    )

    print(
        f"Target-pockets: "
        f"{len(unique_pockets)}"
    )

    print(
        f"Tasks generated: "
        f"{len(jobs)}"
    )

    print(
        f"Job table: {JOB_TABLE}"
    )


if __name__ == "__main__":
    main()
