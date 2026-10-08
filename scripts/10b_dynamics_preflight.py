#!/usr/bin/env python3

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path


ROOT = Path.home() / "PocketToDynamics"

OUT = (
    ROOT /
    "work/dynamics/dynamics_preflight.json"
)

LOG = (
    ROOT /
    "logs/M10B_dynamics_preflight.log"
)


EXPECTED = {
    "gmx":
        "/mnt/labio-fast-disk/"
        "conda-envs/p2d-md/bin/gmx",

    "antechamber":
        str(
            Path.home() /
            "amber/ambertools25/bin/"
            "antechamber"
        ),

    "parmchk2":
        str(
            Path.home() /
            "amber/ambertools25/bin/"
            "parmchk2"
        ),

    "tleap":
        str(
            Path.home() /
            "amber/ambertools25/bin/"
            "tleap"
        ),

    "obabel":
        "/usr/bin/obabel",
}


def version(
    executable: str,
    args: list[str],
):

    try:

        proc = subprocess.run(
            [executable, *args],
            capture_output=True,
            text=True,
            timeout=30,
        )

        text = (
            proc.stdout.strip()
            or
            proc.stderr.strip()
        )

        return {
            "returncode":
                proc.returncode,

            "output":
                text[:4000],
        }

    except Exception as exc:

        return {
            "returncode": -1,
            "output":
                f"{type(exc).__name__}: {exc}",
        }


def main():

    OUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    records = {}

    missing = []

    for name, configured in (
        EXPECTED.items()
    ):

        path = Path(configured)

        found = (
            path.exists()
            and
            path.is_file()
        )

        records[name] = {
            "configured_path":
                configured,

            "found":
                found,
        }

        if not found:
            missing.append(name)

    if records["gmx"]["found"]:

        records["gmx"]["version"] = (
            version(
                EXPECTED["gmx"],
                ["--version"],
            )
        )

    if records[
        "antechamber"
    ]["found"]:

        records[
            "antechamber"
        ]["version"] = version(
            EXPECTED["antechamber"],
            ["-h"],
        )

    if records["tleap"]["found"]:

        records["tleap"]["version"] = (
            version(
                EXPECTED["tleap"],
                ["-h"],
            )
        )

    if records[
        "obabel"
    ]["found"]:

        records[
            "obabel"
        ]["version"] = version(
            EXPECTED["obabel"],
            ["-V"],
        )

    try:

        import parmed

        records["parmed"] = {
            "found": True,
            "version":
                getattr(
                    parmed,
                    "__version__",
                    "unknown",
                ),
        }

    except Exception as exc:

        records["parmed"] = {
            "found": False,
            "error": str(exc),
        }

        missing.append(
            "parmed"
        )

    status = (
        "PASS_M10B_DYNAMICS_PREFLIGHT"
        if not missing
        else
        "FAIL_M10B_DYNAMICS_PREFLIGHT"
    )

    report = {
        "status": status,
        "tools": records,
        "missing": missing,
    }

    OUT.write_text(
        json.dumps(
            report,
            indent=2,
        )
        + "\n"
    )

    LOG.write_text(
        "\n".join(
            [
                status,
                (
                    "Missing: "
                    + (
                        ",".join(missing)
                        if missing
                        else "none"
                    )
                ),
                f"Report: {OUT}",
            ]
        )
        + "\n"
    )

    print(
        "===== M10B DYNAMICS PREFLIGHT ====="
    )

    for name, rec in records.items():

        found = rec.get(
            "found",
            False,
        )

        print(
            f"{name:<14}",
            "FOUND"
            if found
            else "NOT FOUND",
        )

        if rec.get(
            "configured_path"
        ):

            print(
                " " * 15
                + rec[
                    "configured_path"
                ]
            )

    print()
    print(status)

    if missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
