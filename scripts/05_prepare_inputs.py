#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

NORMALIZED = ROOT / "work" / "normalized"

REC_DIR = NORMALIZED / "receptors"
LIG_DIR = NORMALIZED / "ligands"

OUT_ROOT = ROOT / "work" / "prepared"
OUT_REC = OUT_ROOT / "receptors"
OUT_LIG = OUT_ROOT / "ligands"
OUT_LOG = OUT_ROOT / "logs"

MANIFEST_TSV = OUT_ROOT / "preparation_manifest.tsv"
SUMMARY_JSON = OUT_ROOT / "preparation_summary.json"
MODULE_LOG = ROOT / "logs" / "M05_prepare_inputs.log"

MEEKO_REC = shutil.which("mk_prepare_receptor.py")
MEEKO_LIG = shutil.which("mk_prepare_ligand.py")


def sha256(path: Path) -> str:
    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)

    return h.hexdigest()


def run_command(cmd: list[str], stdout_path: Path, stderr_path: Path):
    completed = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    stdout_path.write_text(
        completed.stdout or "",
        errors="replace",
    )

    stderr_path.write_text(
        completed.stderr or "",
        errors="replace",
    )

    return completed


def count_pdbqt_atoms(path: Path) -> int:
    if not path.exists():
        return 0

    count = 0

    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            count += 1

    return count


def prepare_receptor(path: Path) -> dict:
    entity_id = path.stem

    output = OUT_REC / f"{entity_id}.pdbqt"

    stdout = OUT_LOG / f"{entity_id}.receptor.stdout.txt"
    stderr = OUT_LOG / f"{entity_id}.receptor.stderr.txt"

    cmd = [
        MEEKO_REC,
        "--read_pdb",
        str(path),
        "-o",
        str(OUT_REC / entity_id),
        "-p",
        str(output),
        "--compute_charges",
    ]

    completed = run_command(
        cmd,
        stdout,
        stderr,
    )

    atoms = count_pdbqt_atoms(output)

    status = "PASS"

    if completed.returncode != 0:
        status = "FAIL_PREPARATION"

    elif not output.exists():
        status = "FAIL_NO_OUTPUT"

    elif output.stat().st_size == 0:
        status = "FAIL_EMPTY_OUTPUT"

    elif atoms == 0:
        status = "FAIL_NO_ATOMS"

    return {
        "entity_type": "receptor",
        "entity_id": entity_id,
        "input_path": str(path.relative_to(ROOT)),
        "input_sha256": sha256(path),
        "backend": "meeko",
        "command": " ".join(cmd),
        "returncode": completed.returncode,
        "output_path": (
            str(output.relative_to(ROOT))
            if output.exists()
            else ""
        ),
        "output_sha256": (
            sha256(output)
            if output.exists() and output.stat().st_size > 0
            else ""
        ),
        "atom_count": atoms,
        "stdout_path": str(stdout.relative_to(ROOT)),
        "stderr_path": str(stderr.relative_to(ROOT)),
        "status": status,
    }


def prepare_ligand(path: Path) -> dict:
    entity_id = path.stem

    output = OUT_LIG / f"{entity_id}.pdbqt"

    stdout = OUT_LOG / f"{entity_id}.ligand.stdout.txt"
    stderr = OUT_LOG / f"{entity_id}.ligand.stderr.txt"

    cmd = [
        MEEKO_LIG,
        "-i",
        str(path),
        "-o",
        str(output),
    ]

    completed = run_command(
        cmd,
        stdout,
        stderr,
    )

    atoms = count_pdbqt_atoms(output)

    status = "PASS"

    if completed.returncode != 0:
        status = "FAIL_PREPARATION"

    elif not output.exists():
        status = "FAIL_NO_OUTPUT"

    elif output.stat().st_size == 0:
        status = "FAIL_EMPTY_OUTPUT"

    elif atoms == 0:
        status = "FAIL_NO_ATOMS"

    else:
        text = output.read_text(errors="replace")

        if "ROOT" not in text:
            status = "FAIL_MISSING_ROOT"

        elif "TORSDOF" not in text:
            status = "FAIL_MISSING_TORSDOF"

    return {
        "entity_type": "ligand",
        "entity_id": entity_id,
        "input_path": str(path.relative_to(ROOT)),
        "input_sha256": sha256(path),
        "backend": "meeko",
        "command": " ".join(cmd),
        "returncode": completed.returncode,
        "output_path": (
            str(output.relative_to(ROOT))
            if output.exists()
            else ""
        ),
        "output_sha256": (
            sha256(output)
            if output.exists() and output.stat().st_size > 0
            else ""
        ),
        "atom_count": atoms,
        "stdout_path": str(stdout.relative_to(ROOT)),
        "stderr_path": str(stderr.relative_to(ROOT)),
        "status": status,
    }


def main():

    for directory in (
        OUT_ROOT,
        OUT_REC,
        OUT_LIG,
        OUT_LOG,
        MODULE_LOG.parent,
    ):
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    if not MEEKO_REC or not MEEKO_LIG:
        print("FAIL_M05_MEEKO_NOT_FOUND")
        sys.exit(1)

    receptors = sorted(
        p for p in REC_DIR.glob("*")
        if p.is_file()
        and p.suffix.lower() == ".pdb"
    )

    ligand_extensions = {
        ".sdf",
        ".mol2",
        ".mol",
    }

    ligands = sorted(
        p for p in LIG_DIR.glob("*")
        if p.is_file()
        and p.suffix.lower() in ligand_extensions
    )

    rows = []

    for receptor in receptors:
        print(
            f"[M05] receptor {receptor.name}"
        )

        row = prepare_receptor(receptor)
        rows.append(row)

        print(
            f"      {row['status']} "
            f"atoms={row['atom_count']}"
        )

    for ligand in ligands:
        print(
            f"[M05] ligand {ligand.name}"
        )

        row = prepare_ligand(ligand)
        rows.append(row)

        print(
            f"      {row['status']} "
            f"atoms={row['atom_count']}"
        )

    fields = [
        "entity_type",
        "entity_id",
        "input_path",
        "input_sha256",
        "backend",
        "command",
        "returncode",
        "output_path",
        "output_sha256",
        "atom_count",
        "stdout_path",
        "stderr_path",
        "status",
    ]

    with MANIFEST_TSV.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    pass_count = sum(
        row["status"] == "PASS"
        for row in rows
    )

    fail_count = len(rows) - pass_count

    summary = {
        "module": "M05_prepare_inputs",
        "timestamp_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "backend": "meeko",
        "meeko_receptor_cli": MEEKO_REC,
        "meeko_ligand_cli": MEEKO_LIG,
        "receptors": len(receptors),
        "ligands": len(ligands),
        "entities_total": len(rows),
        "pass": pass_count,
        "fail": fail_count,
    }

    SUMMARY_JSON.write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n"
    )

    if not rows:

        status = "SKIP_M05_NO_INPUTS"

    elif fail_count:

        status = "FAIL_M05_INPUT_PREPARATION"

    else:

        status = "PASS_M05_INPUT_PREPARATION"

    MODULE_LOG.write_text(
        "\n".join(
            [
                status,
                f"Receptors: {len(receptors)}",
                f"Ligands: {len(ligands)}",
                f"PASS: {pass_count}",
                f"FAIL: {fail_count}",
                f"Manifest: {MANIFEST_TSV}",
            ]
        )
        + "\n"
    )

    print()
    print(status)
    print(f"Receptors: {len(receptors)}")
    print(f"Ligands: {len(ligands)}")
    print(f"PASS: {pass_count}")
    print(f"FAIL: {fail_count}")
    print(f"Manifest: {MANIFEST_TSV}")

    if fail_count:
        sys.exit(1)


if __name__ == "__main__":
    main()
