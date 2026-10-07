#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import csv
import json
import shutil
import subprocess
import sys

ROOT = Path.home() / "PocketToDynamics"
RECEPTOR_DIR = ROOT / "work/normalized/receptors"
OUT_ROOT = ROOT / "work/pockets"
LOG_DIR = ROOT / "logs"

THRESHOLD = 0.10


def find_concavity():
    return shutil.which("concavity")


def parse_scores(score_file, threshold):
    selected = []

    with score_file.open(errors="replace") as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 3:
                continue

            try:
                residue_id = int(fields[0])
                score = float(fields[2])
            except ValueError:
                continue

            if score > threshold:
                selected.append((residue_id, score))

    if not selected:
        return []

    pockets = []
    current = [selected[0]]

    for residue in selected[1:]:
        previous_residue = current[-1][0]
        current_residue = residue[0]

        if current_residue == previous_residue + 1:
            current.append(residue)
        else:
            pockets.append(current)
            current = [residue]

    pockets.append(current)

    return pockets


def pdb_ca_coordinates(pdb_file):
    coords = {}

    with pdb_file.open(errors="replace") as handle:
        for line in handle:

            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()

            if atom_name != "CA":
                continue

            try:
                residue_number = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except (ValueError, IndexError):
                continue

            coords[residue_number] = (x, y, z)

    return coords


def centroid(residue_ids, coordinate_map):
    xyz = [
        coordinate_map[r]
        for r in residue_ids
        if r in coordinate_map
    ]

    if not xyz:
        return None

    n = len(xyz)

    return (
        sum(p[0] for p in xyz) / n,
        sum(p[1] for p in xyz) / n,
        sum(p[2] for p in xyz) / n,
    )


def process_receptor(receptor, concavity_bin):
    target_id = receptor.stem
    target_dir = OUT_ROOT / target_id
    target_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = target_dir / "concavity"

    cmd = [
        concavity_bin,
        str(receptor),
        str(output_prefix),
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    (target_dir / "concavity_stdout.txt").write_text(
        result.stdout or ""
    )

    (target_dir / "concavity_stderr.txt").write_text(
        result.stderr or ""
    )

    if result.returncode != 0:
        return {
            "target_id": target_id,
            "status": "FAIL_CONCAVITY",
            "pockets": [],
        }

    score_candidates = sorted(
        target_dir.glob("concavity*.scores")
    )

    if not score_candidates:
        return {
            "target_id": target_id,
            "status": "FAIL_NO_SCORE_FILE",
            "pockets": [],
        }

    score_file = score_candidates[0]

    groups = parse_scores(
        score_file,
        THRESHOLD,
    )

    coordinate_map = pdb_ca_coordinates(receptor)

    pocket_records = []

    for idx, group in enumerate(groups, start=1):
        residue_ids = [x[0] for x in group]
        scores = [x[1] for x in group]

        c = centroid(
            residue_ids,
            coordinate_map,
        )

        pocket_records.append({
            "target_id": target_id,
            "pocket_id": f"P{idx:04d}",
            "n_residues": len(residue_ids),
            "first_residue": min(residue_ids),
            "last_residue": max(residue_ids),
            "mean_score": sum(scores) / len(scores),
            "max_score": max(scores),
            "center_x": None if c is None else round(c[0], 3),
            "center_y": None if c is None else round(c[1], 3),
            "center_z": None if c is None else round(c[2], 3),
            "residue_ids": ",".join(
                map(str, residue_ids)
            ),
        })

    output_tsv = target_dir / "pockets.tsv"

    fields = [
        "target_id",
        "pocket_id",
        "n_residues",
        "first_residue",
        "last_residue",
        "mean_score",
        "max_score",
        "center_x",
        "center_y",
        "center_z",
        "residue_ids",
    ]

    with output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(pocket_records)

    summary = {
        "target_id": target_id,
        "status": "PASS",
        "threshold": THRESHOLD,
        "n_pockets": len(pocket_records),
        "score_file": str(
            score_file.relative_to(ROOT)
        ),
        "output": str(
            output_tsv.relative_to(ROOT)
        ),
    }

    with (
        target_dir / "pocket_summary.json"
    ).open("w") as handle:
        json.dump(
            summary,
            handle,
            indent=2,
        )

    return {
        "target_id": target_id,
        "status": "PASS",
        "pockets": pocket_records,
    }


def main():
    OUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    concavity_bin = find_concavity()

    if concavity_bin is None:
        print("FAIL_M03_CONCAVITY_NOT_FOUND")
        sys.exit(2)

    receptors = sorted(
        RECEPTOR_DIR.glob("*.pdb")
    )

    if not receptors:
        log = LOG_DIR / "M03_pocket_mapping.log"

        log.write_text(
            "SKIP_M03_NO_RECEPTORS\n"
            f"Directory: {RECEPTOR_DIR}\n"
        )

        print("SKIP_M03_NO_RECEPTORS")
        print(
            "No normalized PDB receptors found in:",
            RECEPTOR_DIR,
        )

        return

    results = []

    for receptor in receptors:
        print(
            f"[M03] Processing {receptor.name}"
        )

        result = process_receptor(
            receptor,
            concavity_bin,
        )

        results.append(result)

        print(
            f"      status={result['status']} "
            f"pockets={len(result['pockets'])}"
        )

    pass_count = sum(
        r["status"] == "PASS"
        for r in results
    )

    fail_count = (
        len(results) - pass_count
    )

    global_summary = {
        "module": "M03_POCKET_MAPPING",
        "timestamp": datetime.now().isoformat(),
        "threshold": THRESHOLD,
        "targets_total": len(results),
        "targets_pass": pass_count,
        "targets_failed": fail_count,
        "targets": [
            {
                "target_id": r["target_id"],
                "status": r["status"],
                "n_pockets": len(
                    r["pockets"]
                ),
            }
            for r in results
        ],
    }

    summary_path = (
        OUT_ROOT /
        "pocket_mapping_summary.json"
    )

    with summary_path.open("w") as handle:
        json.dump(
            global_summary,
            handle,
            indent=2,
        )

    log_path = (
        LOG_DIR /
        "M03_pocket_mapping.log"
    )

    with log_path.open("w") as handle:
        if fail_count == 0:
            handle.write(
                "PASS_M03_POCKET_MAPPING\n"
            )
        else:
            handle.write(
                "WARN_M03_PARTIAL_FAILURE\n"
            )

        handle.write(
            f"Targets total: {len(results)}\n"
        )
        handle.write(
            f"Targets PASS: {pass_count}\n"
        )
        handle.write(
            f"Targets failed: {fail_count}\n"
        )

    print()
    print("===== M03 SUMMARY =====")
    print(f"Targets: {len(results)}")
    print(f"PASS: {pass_count}")
    print(f"Failed: {fail_count}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
