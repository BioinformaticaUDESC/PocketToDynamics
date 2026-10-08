#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml


ROOT = Path.home() / "PocketToDynamics"

CONFIG = ROOT / "config" / "screening.yaml"

JOB_MATRIX = (
    ROOT /
    "work/jobs/job_matrix.tsv"
)

PREP_MANIFEST = (
    ROOT /
    "work/prepared/preparation_manifest.tsv"
)

SCREENING_ROOT = (
    ROOT /
    "work/screening/raw"
)

MODULE_LOG = (
    ROOT /
    "logs/M05C_screening.log"
)


def sha256(path: Path) -> str:

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def utcnow() -> str:

    return datetime.now(
        timezone.utc
    ).isoformat()


def load_tsv(path: Path):

    with path.open() as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


def load_config():

    with CONFIG.open() as f:

        return yaml.safe_load(f)


def get_prepared_lookup(rows):

    lookup = {}

    for row in rows:

        key = (
            row["entity_type"],
            row["entity_id"],
        )

        lookup[key] = row

    return lookup


def valid_existing_task(
    manifest_path: Path,
    poses_path: Path,
    task_id: str,
    fingerprint: str,
) -> bool:

    if not manifest_path.exists():
        return False

    if not poses_path.exists():
        return False

    if poses_path.stat().st_size == 0:
        return False

    try:

        data = json.loads(
            manifest_path.read_text()
        )

    except Exception:
        return False

    return (
        data.get("fingerprint") == fingerprint
        and data.get("output_sha256") == sha256(poses_path)
        and data.get("task_id") == task_id
        and
        data.get("status") == "PASS"
        and
        data.get("returncode") == 0
    )


def main():

    cfg = load_config()

    vina_cfg = cfg["vina"]
    box_cfg = cfg["box"]
    engine_cfg = cfg["engines"]["vina"]
    exec_cfg = cfg["execution"]

    vina = Path(
        engine_cfg["executable"]
    )

    if not vina.exists():

        print(
            "FAIL_M05C_VINA_NOT_FOUND"
        )

        sys.exit(1)

    version_result = subprocess.run(
        [str(vina), "--version"],
        capture_output=True,
        text=True,
    )

    version_text = (
        version_result.stdout.strip()
        or
        version_result.stderr.strip()
    )

    jobs = load_tsv(
        JOB_MATRIX
    )

    prep_rows = load_tsv(
        PREP_MANIFEST
    )

    prepared = get_prepared_lookup(
        prep_rows
    )

    SCREENING_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    counters = {
        "PASS": 0,
        "FAIL": 0,
        "SKIP_COMPLETED": 0,
    }

    results = []

    for job in jobs:

        task_id = job["task_id"]
        target_id = job["target_id"]
        ligand_id = job["ligand_id"]

        receptor_row = prepared.get(
            ("receptor", target_id)
        )

        ligand_row = prepared.get(
            ("ligand", ligand_id)
        )

        if not receptor_row or not ligand_row:

            counters["FAIL"] += 1

            results.append(
                (
                    task_id,
                    "FAIL_MISSING_PREPARED_INPUT",
                )
            )

            continue

        receptor = (
            ROOT /
            receptor_row["output_path"]
        )

        ligand = (
            ROOT /
            ligand_row["output_path"]
        )

        task_dir = (
            SCREENING_ROOT /
            task_id
        )

        task_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        poses = (
            task_dir /
            "poses.pdbqt"
        )

        vina_log = (
            task_dir /
            "vina.log"
        )

        stdout_path = (
            task_dir /
            "stdout.txt"
        )

        stderr_path = (
            task_dir /
            "stderr.txt"
        )

        manifest_path = (
            task_dir /
            "task_manifest.json"
        )

        fingerprint_data = {
            "job": job, "vina": vina_cfg, "box": box_cfg,
            "engine_path": str(vina.resolve()), "engine_version": version_text,
            "engine_sha256": sha256(vina),
            "receptor_sha256": sha256(receptor), "ligand_sha256": sha256(ligand),
        }
        fingerprint = hashlib.sha256(json.dumps(fingerprint_data, sort_keys=True).encode()).hexdigest()

        if (
            exec_cfg.get("resume", True)
            and
            not exec_cfg.get(
                "overwrite",
                False,
            )
            and
            valid_existing_task(
                manifest_path,
                poses,
                task_id,
                fingerprint,
            )
        ):

            print(
                f"[M05C] {task_id} "
                f"SKIP_COMPLETED"
            )

            counters[
                "SKIP_COMPLETED"
            ] += 1

            results.append(
                (
                    task_id,
                    "SKIP_COMPLETED",
                )
            )

            continue

        # Preserve previous attempts and prevent stale poses from passing a failed rerun.
        old_files = [p for p in [poses, vina_log, stdout_path, stderr_path, manifest_path] if p.exists()]
        if old_files:
            archive = task_dir / "history" / str(time.time_ns())
            archive.mkdir(parents=True)
            for old in old_files:
                shutil.move(str(old), str(archive / old.name))

        command = [
            str(vina),

            "--receptor",
            str(receptor),

            "--ligand",
            str(ligand),

            "--center_x",
            str(job["center_x"]),

            "--center_y",
            str(job["center_y"]),

            "--center_z",
            str(job["center_z"]),

            "--size_x",
            str(box_cfg["size_x"]),

            "--size_y",
            str(box_cfg["size_y"]),

            "--size_z",
            str(box_cfg["size_z"]),

            "--exhaustiveness",
            str(vina_cfg["exhaustiveness"]),

            "--num_modes",
            str(vina_cfg["num_modes"]),

            "--energy_range",
            str(vina_cfg["energy_range"]),

            "--seed",
            str(vina_cfg["seed"]),

            "--out",
            str(poses),
        ]

        print(
            f"[M05C] running {task_id}"
        )

        start_iso = utcnow()
        start = time.time()

        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )

        end = time.time()
        end_iso = utcnow()

        stdout_path.write_text(
            completed.stdout or "",
            errors="replace",
        )

        stderr_path.write_text(
            completed.stderr or "",
            errors="replace",
        )

        # Vina 1.2.x normally reports table to stdout.
        vina_log.write_text(
            completed.stdout or "",
            errors="replace",
        )

        status = "PASS"

        if completed.returncode != 0:

            status = "FAIL_VINA"

        elif not poses.exists():

            status = "FAIL_NO_OUTPUT"

        elif poses.stat().st_size == 0:

            status = "FAIL_EMPTY_OUTPUT"

        manifest = {
            "task_id": task_id,
            "fingerprint": fingerprint,
            "fingerprint_inputs": fingerprint_data,

            "target_id": target_id,
            "pocket_id": job["pocket_id"],
            "ligand_id": ligand_id,

            "engine": "vina",
            "engine_version": version_text,
            "engine_path": str(vina),

            "receptor_path": str(
                receptor.relative_to(ROOT)
            ),

            "ligand_path": str(
                ligand.relative_to(ROOT)
            ),

            "receptor_sha256": (
                sha256(receptor)
                if receptor.exists()
                else None
            ),

            "ligand_sha256": (
                sha256(ligand)
                if ligand.exists()
                else None
            ),

            "box": {
                "center_x": float(
                    job["center_x"]
                ),
                "center_y": float(
                    job["center_y"]
                ),
                "center_z": float(
                    job["center_z"]
                ),
                "size_x": float(
                    box_cfg["size_x"]
                ),
                "size_y": float(
                    box_cfg["size_y"]
                ),
                "size_z": float(
                    box_cfg["size_z"]
                ),
            },

            "vina": {
                "seed": int(
                    vina_cfg["seed"]
                ),
                "exhaustiveness": int(
                    vina_cfg["exhaustiveness"]
                ),
                "num_modes": int(
                    vina_cfg["num_modes"]
                ),
                "energy_range": float(
                    vina_cfg["energy_range"]
                ),
            },

            "command": command,

            "start_time_utc": start_iso,
            "end_time_utc": end_iso,
            "runtime_seconds": round(
                end - start,
                3,
            ),

            "returncode": (
                completed.returncode
            ),

            "output_path": (
                str(
                    poses.relative_to(ROOT)
                )
                if poses.exists()
                else None
            ),

            "output_sha256": (
                sha256(poses)
                if (
                    poses.exists()
                    and
                    poses.stat().st_size > 0
                )
                else None
            ),

            "stdout_path": str(
                stdout_path.relative_to(ROOT)
            ),

            "stderr_path": str(
                stderr_path.relative_to(ROOT)
            ),

            "status": status,
        }

        manifest_path.write_text(
            json.dumps(
                manifest,
                indent=2,
            )
            + "\n"
        )

        if status == "PASS":

            counters["PASS"] += 1

        else:

            counters["FAIL"] += 1

        results.append(
            (
                task_id,
                status,
            )
        )

        print(
            f"        {status} "
            f"rc={completed.returncode}"
        )

    if not jobs:

        module_status = (
            "SKIP_M05C_NO_TASKS"
        )

    elif counters["FAIL"] > 0:

        module_status = (
            "FAIL_M05C_SCREENING"
        )

    else:

        module_status = (
            "PASS_M05C_SCREENING"
        )

    MODULE_LOG.write_text(
        "\n".join(
            [
                module_status,
                f"Tasks total: {len(jobs)}",
                f"PASS: {counters['PASS']}",
                (
                    "SKIP_COMPLETED: "
                    f"{counters['SKIP_COMPLETED']}"
                ),
                f"FAIL: {counters['FAIL']}",
                f"Engine: {version_text}",
            ]
        )
        + "\n"
    )

    print()
    print(module_status)
    print(
        f"Tasks total: {len(jobs)}"
    )
    print(
        f"PASS: {counters['PASS']}"
    )
    print(
        "SKIP_COMPLETED:",
        counters["SKIP_COMPLETED"],
    )
    print(
        f"FAIL: {counters['FAIL']}"
    )

    if counters["FAIL"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
