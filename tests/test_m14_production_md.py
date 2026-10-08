#!/usr/bin/env python3

import json
from pathlib import Path


JOB = Path(
    "/mnt/labio-fast-disk/"
    "PocketToDynamics/dynamics/"
    "production_fixture/"
    "MDTEST_a5b7f2f48479"
)

manifest = (
    JOB /
    "production_manifest.json"
)

assert manifest.exists(), manifest

data = json.loads(
    manifest.read_text()
)

assert data["status"] == "PASS"
assert data["profile"] == "fixture"

assert data["target_duration_ps"] == 20.0
assert data["completed_duration_ps"] == 20.0

assert data["segments_total"] == 2
assert data["segments_completed"] == 2

assert data["timestep_ps"] == 0.002

segments = data["segments"]

assert len(segments) == 2

assert segments[0]["start_ps"] == 0.0
assert segments[0]["end_ps"] == 10.0

assert segments[1]["start_ps"] == 10.0
assert segments[1]["end_ps"] == 20.0


for index in [1, 2]:

    seg = (
        JOB /
        f"segment_{index:04d}"
    )

    marker = (
        seg /
        "segment_complete.json"
    )

    assert marker.exists(), marker

    marker_data = json.loads(
        marker.read_text()
    )

    assert marker_data["status"] == "PASS"

    for name in [
        "md.gro",
        "md.cpt",
        "md.edr",
        "md.log",
        "md.xtc",
        "md.tpr",
    ]:

        path = seg / name

        assert path.exists(), path
        assert path.stat().st_size > 0


warning_patterns = [
    "LINCS WARNING",
    "constraint warning",
    "nan occurred",
    "fatal error",
    "segmentation fault",
    "blowing up",
]

for index in [1, 2]:

    seg = (
        JOB /
        f"segment_{index:04d}"
    )

    for filename in [
        "md.log",
        "mdrun.stderr",
    ]:

        text = (
            seg /
            filename
        ).read_text(
            errors="replace"
        ).lower()

        for pattern in warning_patterns:

            assert (
                pattern.lower()
                not in text
            ), (
                f"{pattern} found "
                f"in segment "
                f"{index}: "
                f"{filename}"
            )


print(
    "PASS_M14_PRODUCTION_MD_QC"
)

print(
    "Completed duration:",
    data["completed_duration_ps"],
    "ps",
)

print(
    "Segments:",
    data["segments_completed"],
)

print(
    "Final GRO:",
    data["files"]["final_gro"],
)

print(
    "Final CPT:",
    data["files"]["final_cpt"],
)

# Exercise fixture checkpoint behavior before downstream analysis consumes the manifest.
import subprocess
import sys
subprocess.run([sys.executable, str(Path(__file__).with_name('test_m14_resume_integration.py'))], check=True)
