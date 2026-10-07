#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import hashlib
import csv
import json
import re
import shutil
import sys

ROOT = Path.home() / "PocketToDynamics"

INPUT_RECEPTORS = ROOT / "input/receptors"
INPUT_LIGANDS = ROOT / "input/ligands"

OUT_ROOT = ROOT / "work/normalized"
OUT_RECEPTORS = OUT_ROOT / "receptors"
OUT_LIGANDS = OUT_ROOT / "ligands"

LOG_DIR = ROOT / "logs"
MANIFEST_DIR = ROOT / "manifest"

SUPPORTED_RECEPTORS = {".pdb", ".cif", ".mmcif"}
SUPPORTED_LIGANDS = {".sdf", ".mol2", ".pdb", ".smi", ".smiles"}

def sha256_file(path, block_size=1024 * 1024):
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(block_size):
            h.update(chunk)
    return h.hexdigest()

def safe_id(name):
    stem = Path(name).stem

    stem = stem.strip()
    stem = re.sub(r"\s+", "_", stem)
    stem = re.sub(r"[^A-Za-z0-9._-]", "", stem)
    stem = re.sub(r"_+", "_", stem)

    return stem.strip("_")

def scan_inputs(folder, valid_ext):
    files = []

    for path in sorted(folder.rglob("*")):
        if not path.is_file():
            continue

        ext = path.suffix.lower()

        if ext not in valid_ext:
            continue

        files.append(path)

    return files

def normalize_group(files, destination, prefix):
    destination.mkdir(parents=True, exist_ok=True)

    records = []
    used_ids = {}

    for idx, src in enumerate(files, start=1):

        base = safe_id(src.name)

        if not base:
            base = f"{prefix}_{idx:05d}"

        internal_id = base

        if internal_id in used_ids:
            used_ids[internal_id] += 1
            internal_id = f"{internal_id}_{used_ids[internal_id]}"
        else:
            used_ids[internal_id] = 1

        ext = src.suffix.lower()

        dst_name = f"{internal_id}{ext}"
        dst = destination / dst_name

        shutil.copy2(src, dst)

        records.append({
            "internal_id": internal_id,
            "original_filename": src.name,
            "original_path": str(src.relative_to(ROOT)),
            "normalized_path": str(dst.relative_to(ROOT)),
            "extension": ext,
            "size_bytes": dst.stat().st_size,
            "sha256": sha256_file(dst),
        })

    return records

def write_tsv(records, path):
    fields = [
        "internal_id",
        "original_filename",
        "original_path",
        "normalized_path",
        "extension",
        "size_bytes",
        "sha256",
    ]

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for row in records:
            writer.writerow(row)

def main():

    OUT_RECEPTORS.mkdir(parents=True, exist_ok=True)
    OUT_LIGANDS.mkdir(parents=True, exist_ok=True)

    receptor_files = scan_inputs(
        INPUT_RECEPTORS,
        SUPPORTED_RECEPTORS,
    )

    ligand_files = scan_inputs(
        INPUT_LIGANDS,
        SUPPORTED_LIGANDS,
    )

    receptor_records = normalize_group(
        receptor_files,
        OUT_RECEPTORS,
        "REC",
    )

    ligand_records = normalize_group(
        ligand_files,
        OUT_LIGANDS,
        "LIG",
    )

    write_tsv(
        receptor_records,
        OUT_ROOT / "receptor_id_map.tsv",
    )

    write_tsv(
        ligand_records,
        OUT_ROOT / "ligand_id_map.tsv",
    )

    summary = {
        "module": "M02_INPUT_NORMALIZATION",
        "timestamp": datetime.now().isoformat(),
        "receptors_detected": len(receptor_records),
        "ligands_detected": len(ligand_records),
        "receptor_extensions": sorted(
            set(x["extension"] for x in receptor_records)
        ),
        "ligand_extensions": sorted(
            set(x["extension"] for x in ligand_records)
        ),
    }

    summary_path = OUT_ROOT / "normalization_summary.json"

    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)

    log_path = LOG_DIR / "M02_normalization.log"

    with log_path.open("w") as f:

        if not receptor_records:
            f.write("WARN_NO_RECEPTORS\n")

        if not ligand_records:
            f.write("WARN_NO_LIGANDS\n")

        f.write("PASS_M02_INPUT_NORMALIZATION\n")
        f.write(
            f"Receptors: {len(receptor_records)}\n"
        )
        f.write(
            f"Ligands: {len(ligand_records)}\n"
        )

    print("PASS_M02_INPUT_NORMALIZATION")
    print(f"Receptors detected: {len(receptor_records)}")
    print(f"Ligands detected: {len(ligand_records)}")
    print(f"Output: {OUT_ROOT}")

if __name__ == "__main__":
    main()

