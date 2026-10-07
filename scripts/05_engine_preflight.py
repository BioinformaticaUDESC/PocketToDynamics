#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import json
import shutil
import subprocess

ROOT = Path.home() / "PocketToDynamics"

OUT_DIR = ROOT / "work/screening"
LOG_DIR = ROOT / "logs"

SUMMARY_FILE = OUT_DIR / "engine_preflight.json"
LOG_FILE = LOG_DIR / "M05_engine_preflight.log"


TOOLS = {
    "vina": [
        "vina",
    ],

    "gnina": [
        "gnina",
    ],

    "smina": [
        "smina",
    ],

    "obabel": [
        "obabel",
        "babel",
    ],

    "pythonsh": [
        "pythonsh",
    ],

    "prepare_receptor4": [
        "prepare_receptor4.py",
    ],

    "prepare_ligand4": [
        "prepare_ligand4.py",
    ],
}


def first_executable(names):
    for name in names:
        path = shutil.which(name)

        if path:
            return path

    return None


def get_version(path):
    commands = [
        [path, "--version"],
        [path, "-V"],
        [path, "-v"],
    ]

    for cmd in commands:
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=10,
            )

        except Exception:
            continue

        text = (
            (result.stdout or "")
            + "\n"
            + (result.stderr or "")
        ).strip()

        if text:
            return text.splitlines()[0][:300]

    return None


def main():
    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    LOG_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    detected = {}

    for tool, candidates in TOOLS.items():
        path = first_executable(
            candidates
        )

        detected[tool] = {
            "available": path is not None,
            "path": path,
            "version": (
                get_version(path)
                if path
                else None
            ),
        }

    engines = [
        name
        for name in (
            "vina",
            "gnina",
            "smina",
        )
        if detected[name]["available"]
    ]

    preparation_tools = {
        "openbabel":
            detected["obabel"]["available"],

        "mgltools_receptor":
            detected[
                "prepare_receptor4"
            ]["available"],

        "mgltools_ligand":
            detected[
                "prepare_ligand4"
            ]["available"],
    }

    if engines:
        status = (
            "PASS_M05_ENGINE_PREFLIGHT"
        )
    else:
        status = (
            "WARN_M05_NO_SCREENING_ENGINE"
        )

    summary = {
        "module":
            "M05_ENGINE_PREFLIGHT",

        "timestamp":
            datetime.now().isoformat(),

        "status":
            status,

        "screening_engines":
            engines,

        "n_screening_engines":
            len(engines),

        "preparation_tools":
            preparation_tools,

        "tools":
            detected,
    }

    with SUMMARY_FILE.open(
        "w"
    ) as handle:
        json.dump(
            summary,
            handle,
            indent=2,
        )

    with LOG_FILE.open(
        "w"
    ) as handle:
        handle.write(
            status + "\n"
        )

        handle.write(
            "Screening engines: "
            + (
                ",".join(engines)
                if engines
                else "none"
            )
            + "\n"
        )

        for tool, data in detected.items():
            handle.write(
                f"{tool}: "
                f"{'FOUND' if data['available'] else 'NOT_FOUND'}"
            )

            if data["path"]:
                handle.write(
                    f" | {data['path']}"
                )

            if data["version"]:
                handle.write(
                    f" | {data['version']}"
                )

            handle.write("\n")

    print(
        "===== M05 ENGINE PREFLIGHT ====="
    )

    print(
        f"Status: {status}"
    )

    print()

    for tool, data in detected.items():
        flag = (
            "FOUND"
            if data["available"]
            else "NOT FOUND"
        )

        print(
            f"{tool:20s} {flag}"
        )

        if data["path"]:
            print(
                f"{'':20s} {data['path']}"
            )

        if data["version"]:
            print(
                f"{'':20s} {data['version']}"
            )

    print()
    print(
        "Screening engines:",
        ", ".join(engines)
        if engines
        else "none",
    )

    print(
        f"Summary: {SUMMARY_FILE}"
    )


if __name__ == "__main__":
    main()
