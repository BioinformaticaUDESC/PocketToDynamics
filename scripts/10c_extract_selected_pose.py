#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

DEFAULT_MANIFEST = (
    ROOT /
    "work/dynamics/jobs/"
    "dynamics_job_manifest.tsv"
)

DEFAULT_OUTPUT = (
    ROOT /
    "work/dynamics/prepared"
)

DEFAULT_LOG = (
    ROOT /
    "logs/M10C_pose_extraction.log"
)

OBABEL = Path("/usr/bin/obabel")


def sha256(path: Path) -> str:

    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


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


def extract_model(
    source: Path,
    mode: int,
) -> str:

    lines = source.read_text(
        errors="replace"
    ).splitlines()

    selected = []
    inside = False
    saw_models = False

    for line in lines:

        stripped = line.strip()

        if stripped.startswith("MODEL"):

            saw_models = True

            parts = stripped.split()

            try:
                model_number = int(
                    parts[1]
                )
            except Exception:
                continue

            if model_number == mode:
                inside = True
                selected.append(line)

            else:
                inside = False

            continue

        if inside:

            selected.append(line)

            if stripped == "ENDMDL":
                break

    if selected:
        return "\n".join(selected) + "\n"

    # Defensive fallback for single-pose PDBQT
    if not saw_models and mode == 1:

        return (
            "\n".join(lines)
            + "\n"
        )

    raise ValueError(
        f"MODEL {mode} not found in {source}"
    )


def run_obabel(
    source: Path,
    output: Path,
    fmt: str,
):

    command = [
        str(OBABEL),
        "-ipdbqt",
        str(source),
        f"-o{fmt}",
        "-O",
        str(output),
    ]

    proc = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    return (
        command,
        proc,
    )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT,
    )

    parser.add_argument(
        "--log",
        type=Path,
        default=DEFAULT_LOG,
    )

    args = parser.parse_args()

    manifest_path = (
        args.manifest.expanduser().resolve()
    )

    output_root = (
        args.output_root.expanduser().resolve()
    )

    log_path = (
        args.log.expanduser().resolve()
    )

    output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    log_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    jobs = load_tsv(
        manifest_path
    )

    summary_table = (
        output_root /
        "pose_extraction_manifest.tsv"
    )

    fields = [
        "dynamics_job_id",
        "target_id",
        "ligand_id",
        "pocket_id",
        "screening_task_id",
        "selected_mode",
        "selected_score",
        "selected_pose_pdbqt",
        "selected_pose_mol2",
        "selected_pose_pdb",
        "receptor_pdb",
        "status",
    ]

    if not jobs:

        with summary_table.open(
            "w",
            newline="",
        ) as f:

            writer = csv.DictWriter(
                f,
                fieldnames=fields,
                delimiter="\t",
            )

            writer.writeheader()

        status = (
            "SKIP_M10C_NO_DYNAMICS_JOBS"
        )

        log_path.write_text(
            "\n".join(
                [
                    status,
                    "Dynamics jobs: 0",
                    "Extracted poses: 0",
                ]
            )
            + "\n"
        )

        print(status)
        print("Dynamics jobs: 0")
        print("Extracted poses: 0")

        return

    if not OBABEL.exists():

        print(
            "FAIL_M10C_OBABEL_NOT_FOUND"
        )

        raise SystemExit(1)

    results = []
    failures = 0

    for job in jobs:

        job_id = (
            job["dynamics_job_id"]
        )

        target_id = (
            job["target_id"]
        )

        ligand_id = (
            job["ligand_id"]
        )

        mode = int(
            job["selected_mode"]
        )

        source = (
            ROOT /
            job["source_pose_file"]
        )

        receptor = (
            ROOT /
            "work/normalized/receptors" /
            f"{target_id}.pdb"
        )

        job_dir = (
            output_root /
            job_id
        )

        job_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        selected_pdbqt = (
            job_dir /
            "selected_pose.pdbqt"
        )

        selected_mol2 = (
            job_dir /
            "selected_pose.mol2"
        )

        selected_pdb = (
            job_dir /
            "selected_pose.pdb"
        )

        receptor_copy = (
            job_dir /
            "receptor.pdb"
        )

        pose_manifest = (
            job_dir /
            "pose_manifest.json"
        )

        try:

            if not source.exists():
                raise FileNotFoundError(
                    f"Missing source pose: {source}"
                )

            if not receptor.exists():
                raise FileNotFoundError(
                    f"Missing receptor: {receptor}"
                )

            selected_text = (
                extract_model(
                    source,
                    mode,
                )
            )

            selected_pdbqt.write_text(
                selected_text
            )

            command_mol2, p_mol2 = (
                run_obabel(
                    selected_pdbqt,
                    selected_mol2,
                    "mol2",
                )
            )

            if (
                p_mol2.returncode != 0
                or
                not selected_mol2.exists()
                or
                selected_mol2.stat().st_size == 0
            ):
                raise RuntimeError(
                    "OpenBabel MOL2 conversion failed: "
                    + p_mol2.stderr
                )

            command_pdb, p_pdb = (
                run_obabel(
                    selected_pdbqt,
                    selected_pdb,
                    "pdb",
                )
            )

            if (
                p_pdb.returncode != 0
                or
                not selected_pdb.exists()
                or
                selected_pdb.stat().st_size == 0
            ):
                raise RuntimeError(
                    "OpenBabel PDB conversion failed: "
                    + p_pdb.stderr
                )

            shutil.copy2(
                receptor,
                receptor_copy,
            )

            metadata = {
                "dynamics_job_id":
                    job_id,

                "target_id":
                    target_id,

                "ligand_id":
                    ligand_id,

                "pocket_id":
                    job["pocket_id"],

                "screening_task_id":
                    job[
                        "screening_task_id"
                    ],

                "selected_mode":
                    mode,

                "selected_score":
                    float(
                        job["selected_score"]
                    ),

                "source_pose_file":
                    str(
                        source.relative_to(
                            ROOT
                        )
                    ),

                "source_pose_sha256":
                    sha256(source),

                "selected_pose_pdbqt":
                    str(
                        selected_pdbqt.relative_to(
                            ROOT
                        )
                    )
                    if selected_pdbqt.is_relative_to(ROOT)
                    else str(selected_pdbqt),

                "selected_pose_sha256":
                    sha256(
                        selected_pdbqt
                    ),

                "selected_mol2_sha256":
                    sha256(
                        selected_mol2
                    ),

                "selected_pdb_sha256":
                    sha256(
                        selected_pdb
                    ),

                "receptor_sha256":
                    sha256(
                        receptor_copy
                    ),

                "obabel_mol2_command":
                    command_mol2,

                "obabel_pdb_command":
                    command_pdb,

                "status":
                    "PASS",
            }

            pose_manifest.write_text(
                json.dumps(
                    metadata,
                    indent=2,
                )
                + "\n"
            )

            result = {
                "dynamics_job_id":
                    job_id,

                "target_id":
                    target_id,

                "ligand_id":
                    ligand_id,

                "pocket_id":
                    job["pocket_id"],

                "screening_task_id":
                    job[
                        "screening_task_id"
                    ],

                "selected_mode":
                    mode,

                "selected_score":
                    job[
                        "selected_score"
                    ],

                "selected_pose_pdbqt":
                    str(selected_pdbqt),

                "selected_pose_mol2":
                    str(selected_mol2),

                "selected_pose_pdb":
                    str(selected_pdb),

                "receptor_pdb":
                    str(receptor_copy),

                "status":
                    "PASS",
            }

            results.append(
                result
            )

            print(
                f"[M10C] {job_id} "
                f"PASS mode={mode}"
            )

        except Exception as exc:

            failures += 1

            results.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "target_id":
                        target_id,

                    "ligand_id":
                        ligand_id,

                    "pocket_id":
                        job["pocket_id"],

                    "screening_task_id":
                        job[
                            "screening_task_id"
                        ],

                    "selected_mode":
                        mode,

                    "selected_score":
                        job[
                            "selected_score"
                        ],

                    "selected_pose_pdbqt":
                        "",

                    "selected_pose_mol2":
                        "",

                    "selected_pose_pdb":
                        "",

                    "receptor_pdb":
                        "",

                    "status":
                        "FAIL",
                }
            )

            print(
                f"[M10C] {job_id} "
                f"FAIL {type(exc).__name__}: "
                f"{exc}"
            )

    with summary_table.open(
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
            results
        )

    passed = sum(
        row["status"] == "PASS"
        for row in results
    )

    if failures:

        status = (
            "FAIL_M10C_POSE_EXTRACTION"
        )

    else:

        status = (
            "PASS_M10C_POSE_EXTRACTION"
        )

    log_path.write_text(
        "\n".join(
            [
                status,
                f"Dynamics jobs: {len(jobs)}",
                f"Extracted poses: {passed}",
                f"Failures: {failures}",
                (
                    "Manifest: "
                    f"{summary_table}"
                ),
            ]
        )
        + "\n"
    )

    print()
    print(status)
    print(
        "Dynamics jobs:",
        len(jobs),
    )
    print(
        "Extracted poses:",
        passed,
    )
    print(
        "Failures:",
        failures,
    )

    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
