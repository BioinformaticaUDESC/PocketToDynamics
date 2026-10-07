#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import csv
import re

ROOT = Path.home() / "PocketToDynamics"
LOG_DIR = ROOT / "logs"
STATUS_FILE = ROOT / "manifest" / "STATUS.tsv"

MODULE_LOGS = [
    ("M01", "M01_init.log"),
    ("M02", "M02_normalization.log"),
    ("M03", "M03_pocket_mapping.log"),
    ("M04", "M04_job_builder.log"),
]

def classify(line):
    line = line.strip().upper()

    if line.startswith("PASS"):
        return "PASS"

    if line.startswith("WARN"):
        return "WARN"

    if line.startswith("FAIL") or line.startswith("ERROR"):
        return "FAIL"

    if line.startswith("SKIP"):
        return "SKIP"

    return "UNKNOWN"

def main():
    rows = []

    for module, filename in MODULE_LOGS:
        path = LOG_DIR / filename

        if not path.exists():
            rows.append({
                "module": module,
                "status": "NOT_RUN",
                "message": "",
                "log_file": str(path.relative_to(ROOT)),
            })
            continue

        lines = path.read_text(errors="replace").splitlines()

        first = lines[0] if lines else ""

        rows.append({
            "module": module,
            "status": classify(first),
            "message": first,
            "log_file": str(path.relative_to(ROOT)),
        })

    with STATUS_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "module",
                "status",
                "message",
                "log_file",
            ],
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    counts = {}

    for row in rows:
        counts[row["status"]] = counts.get(
            row["status"],
            0
        ) + 1

    print("===== PIPELINE STATUS =====")

    for row in rows:
        print(
            f"{row['module']:4s} "
            f"{row['status']:8s} "
            f"{row['message']}"
        )

    print()
    print("Summary:", counts)
    print(f"Status file: {STATUS_FILE}")

if __name__ == "__main__":
    main()
