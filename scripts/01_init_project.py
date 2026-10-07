#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import hashlib
import json
import platform
import sys

ROOT = Path.home() / "PocketToDynamics"

REQUIRED_DIRS = [
    "config",
    "input/receptors",
    "input/ligands",
    "input/metadata",
    "scripts",
    "work",
    "results/tables",
    "results/figures",
    "results/reports",
    "results/supplementary",
    "logs",
    "manifest",
]

def sha256_file(path, block_size=1024 * 1024):
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(block_size):
            h.update(chunk)
    return h.hexdigest()

def main():
    for d in REQUIRED_DIRS:
        (ROOT / d).mkdir(parents=True, exist_ok=True)

    inputs = []

    for subdir in [
        ROOT / "input/receptors",
        ROOT / "input/ligands",
        ROOT / "input/metadata",
    ]:
        for f in sorted(subdir.rglob("*")):
            if f.is_file():
                inputs.append({
                    "path": str(f.relative_to(ROOT)),
                    "size_bytes": f.stat().st_size,
                    "sha256": sha256_file(f),
                })

    manifest = {
        "project": "PocketToDynamics",
        "created_at": datetime.now().isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "root": str(ROOT),
        "n_input_files": len(inputs),
        "inputs": inputs,
    }

    out = ROOT / "manifest/project_manifest.json"

    with out.open("w") as f:
        json.dump(manifest, f, indent=2)

    log = ROOT / "logs/M01_init.log"

    with log.open("w") as f:
        f.write("PASS_M01_PROJECT_INIT\n")
        f.write(f"Root: {ROOT}\n")
        f.write(f"Input files: {len(inputs)}\n")
        f.write(f"Manifest: {out}\n")

    print("PASS_M01_PROJECT_INIT")
    print(f"Project root: {ROOT}")
    print(f"Input files detected: {len(inputs)}")
    print(f"Manifest: {out}")

if __name__ == "__main__":
    main()
