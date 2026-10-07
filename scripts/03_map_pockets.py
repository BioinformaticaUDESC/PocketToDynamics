#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
from math import sqrt
import csv
import json
import shutil
import subprocess
import sys

ROOT = Path.home() / "PocketToDynamics"

RECEPTOR_DIR = ROOT / "work/normalized/receptors"
OUT_ROOT = ROOT / "work/pockets"
LOG_DIR = ROOT / "logs"

SCORE_THRESHOLD = 0.10
SPATIAL_CUTOFF = 8.0
MIN_RESIDUES = 3


def find_concavity():
    return shutil.which("concavity")


def distance(a, b):
    return sqrt(
        (a[0] - b[0]) ** 2 +
        (a[1] - b[1]) ** 2 +
        (a[2] - b[2]) ** 2
    )


def parse_pdb_ca(pdb_file):
    residues = {}

    with pdb_file.open(errors="replace") as handle:
        for line in handle:

            if not line.startswith("ATOM"):
                continue

            if line[12:16].strip() != "CA":
                continue

            chain = line[21].strip() or "_"

            try:
                resnum = int(
                    line[22:26].strip()
                )

                x = float(
                    line[30:38].strip()
                )

                y = float(
                    line[38:46].strip()
                )

                z = float(
                    line[46:54].strip()
                )

            except ValueError:
                continue

            residues[
                (chain, resnum)
            ] = (x, y, z)

    return residues


def chain_from_score_filename(path, target_id):
    name = path.name

    prefix = f"{target_id}_"
    suffix = "_concavity.scores"

    if (
        name.startswith(prefix)
        and name.endswith(suffix)
    ):
        chain = name[
            len(prefix):
            -len(suffix)
        ]

        return chain or "_"

    return None


def parse_concavity_scores(
    score_files,
    target_id,
    coordinates,
):
    candidates = []

    for score_file in score_files:

        chain = chain_from_score_filename(
            score_file,
            target_id,
        )

        if chain is None:
            continue

        with score_file.open(
            errors="replace"
        ) as handle:

            for line in handle:

                line = line.strip()

                if (
                    not line
                    or line.startswith("#")
                ):
                    continue

                fields = line.split()

                if len(fields) < 3:
                    continue

                try:
                    resnum = int(fields[0])
                    score = float(fields[2])
                except ValueError:
                    continue

                key = (
                    chain,
                    resnum,
                )

                if (
                    score >= SCORE_THRESHOLD
                    and key in coordinates
                ):
                    candidates.append({
                        "chain":
                            chain,

                        "resnum":
                            resnum,

                        "score":
                            score,

                        "coord":
                            coordinates[key],
                    })

    return candidates


def spatial_components(candidates):
    """
    Connected-component clustering in Cartesian space.

    Residues are connected when their C-alpha
    atoms are within SPATIAL_CUTOFF Angstrom.
    """

    n = len(candidates)

    adjacency = [
        []
        for _ in range(n)
    ]

    for i in range(n):

        for j in range(i + 1, n):

            if distance(
                candidates[i]["coord"],
                candidates[j]["coord"],
            ) <= SPATIAL_CUTOFF:

                adjacency[i].append(j)
                adjacency[j].append(i)

    visited = set()
    components = []

    for start in range(n):

        if start in visited:
            continue

        stack = [start]
        visited.add(start)

        component = []

        while stack:

            node = stack.pop()

            component.append(
                candidates[node]
            )

            for neighbour in adjacency[node]:

                if neighbour not in visited:

                    visited.add(
                        neighbour
                    )

                    stack.append(
                        neighbour
                    )

        components.append(
            component
        )

    return components


def centroid(component):
    xyz = [
        r["coord"]
        for r in component
    ]

    n = len(xyz)

    return (
        sum(x[0] for x in xyz) / n,
        sum(x[1] for x in xyz) / n,
        sum(x[2] for x in xyz) / n,
    )


def process_receptor(
    receptor,
    concavity_bin,
):
    target_id = receptor.stem

    target_dir = (
        OUT_ROOT /
        target_id
    )

    target_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    #
    # Remove only generated Concavity/runtime files.
    #
    for old in target_dir.iterdir():

        if old.is_file():
            old.unlink()

    cmd = [
        concavity_bin,
        str(receptor.resolve()),
        "concavity",
    ]

    result = subprocess.run(
        cmd,
        cwd=target_dir,
        capture_output=True,
        text=True,
    )

    (
        target_dir /
        "concavity_stdout.txt"
    ).write_text(
        result.stdout or ""
    )

    (
        target_dir /
        "concavity_stderr.txt"
    ).write_text(
        result.stderr or ""
    )

    if result.returncode != 0:

        return {
            "target_id":
                target_id,

            "status":
                "FAIL_CONCAVITY",

            "pockets":
                [],
        }

    score_files = sorted(
        target_dir.glob(
            "*_concavity.scores"
        )
    )

    if not score_files:

        return {
            "target_id":
                target_id,

            "status":
                "FAIL_NO_SCORE_FILES",

            "pockets":
                [],
        }

    coordinates = parse_pdb_ca(
        receptor
    )

    candidates = (
        parse_concavity_scores(
            score_files,
            target_id,
            coordinates,
        )
    )

    if not candidates:

        return {
            "target_id":
                target_id,

            "status":
                "WARN_NO_CANDIDATE_RESIDUES",

            "pockets":
                [],
        }

    components = spatial_components(
        candidates
    )

    #
    # Remove isolated / very small components.
    #
    components = [
        c
        for c in components
        if len(c) >= MIN_RESIDUES
    ]

    #
    # Strongest / largest pockets first.
    #
    components.sort(
        key=lambda c: (
            -sum(
                r["score"]
                for r in c
            ),
            -len(c),
        )
    )

    pocket_records = []

    for index, component in enumerate(
        components,
        start=1,
    ):

        c = centroid(
            component
        )

        scores = [
            r["score"]
            for r in component
        ]

        chains = sorted({
            r["chain"]
            for r in component
        })

        residues = sorted(
            component,
            key=lambda r: (
                r["chain"],
                r["resnum"],
            ),
        )

        residue_labels = [
            f"{r['chain']}:{r['resnum']}"
            for r in residues
        ]

        pocket_records.append({
            "target_id":
                target_id,

            "pocket_id":
                f"P{index:04d}",

            "n_residues":
                len(component),

            "chains":
                ",".join(chains),

            "mean_score":
                round(
                    sum(scores) /
                    len(scores),
                    6,
                ),

            "max_score":
                round(
                    max(scores),
                    6,
                ),

            "sum_score":
                round(
                    sum(scores),
                    6,
                ),

            "center_x":
                round(c[0], 3),

            "center_y":
                round(c[1], 3),

            "center_z":
                round(c[2], 3),

            "residues":
                ",".join(
                    residue_labels
                ),
        })

    output_tsv = (
        target_dir /
        "pockets.tsv"
    )

    fields = [
        "target_id",
        "pocket_id",
        "n_residues",
        "chains",
        "mean_score",
        "max_score",
        "sum_score",
        "center_x",
        "center_y",
        "center_z",
        "residues",
    ]

    with output_tsv.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        writer.writerows(
            pocket_records
        )

    if pocket_records:
        status = "PASS"
    else:
        status = "WARN_NO_SPATIAL_POCKETS"

    summary = {
        "target_id":
            target_id,

        "status":
            status,

        "score_threshold":
            SCORE_THRESHOLD,

        "spatial_cutoff_angstrom":
            SPATIAL_CUTOFF,

        "minimum_residues":
            MIN_RESIDUES,

        "score_files": [
            str(
                f.relative_to(ROOT)
            )
            for f in score_files
        ],

        "n_score_files":
            len(score_files),

        "n_candidate_residues":
            len(candidates),

        "n_spatial_components_raw":
            len(
                spatial_components(
                    candidates
                )
            ),

        "n_pockets":
            len(pocket_records),

        "pocket_table":
            str(
                output_tsv.relative_to(
                    ROOT
                )
            ),
    }

    with (
        target_dir /
        "pocket_summary.json"
    ).open("w") as handle:

        json.dump(
            summary,
            handle,
            indent=2,
        )

    return {
        "target_id":
            target_id,

        "status":
            status,

        "pockets":
            pocket_records,
    }


def main():

    OUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    concavity_bin = (
        find_concavity()
    )

    if concavity_bin is None:

        print(
            "FAIL_M03_CONCAVITY_NOT_FOUND"
        )

        sys.exit(2)

    receptors = sorted(
        RECEPTOR_DIR.glob(
            "*.pdb"
        )
    )

    if not receptors:

        log = (
            LOG_DIR /
            "M03_pocket_mapping.log"
        )

        log.write_text(
            "SKIP_M03_NO_RECEPTORS\n"
            f"Directory: "
            f"{RECEPTOR_DIR}\n"
        )

        print(
            "SKIP_M03_NO_RECEPTORS"
        )

        return

    results = []

    for receptor in receptors:

        print(
            f"[M03] Processing "
            f"{receptor.name}"
        )

        result = process_receptor(
            receptor,
            concavity_bin,
        )

        results.append(
            result
        )

        print(
            f"      status="
            f"{result['status']} "
            f"pockets="
            f"{len(result['pockets'])}"
        )

    pass_count = sum(
        r["status"] == "PASS"
        for r in results
    )

    warn_count = sum(
        r["status"].startswith(
            "WARN"
        )
        for r in results
    )

    fail_count = sum(
        r["status"].startswith(
            "FAIL"
        )
        for r in results
    )

    summary = {
        "module":
            "M03_POCKET_MAPPING",

        "timestamp":
            datetime.now().isoformat(),

        "method":
            "Concavity score threshold + "
            "C-alpha spatial clustering",

        "score_threshold":
            SCORE_THRESHOLD,

        "spatial_cutoff_angstrom":
            SPATIAL_CUTOFF,

        "minimum_residues":
            MIN_RESIDUES,

        "targets_total":
            len(results),

        "targets_pass":
            pass_count,

        "targets_warn":
            warn_count,

        "targets_failed":
            fail_count,

        "targets": [
            {
                "target_id":
                    r["target_id"],

                "status":
                    r["status"],

                "n_pockets":
                    len(
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

    with summary_path.open(
        "w"
    ) as handle:

        json.dump(
            summary,
            handle,
            indent=2,
        )

    if fail_count:
        module_status = (
            "WARN_M03_PARTIAL_FAILURE"
        )

    elif warn_count:
        module_status = (
            "WARN_M03_POCKET_WARNINGS"
        )

    else:
        module_status = (
            "PASS_M03_POCKET_MAPPING"
        )

    log_path = (
        LOG_DIR /
        "M03_pocket_mapping.log"
    )

    with log_path.open(
        "w"
    ) as handle:

        handle.write(
            module_status +
            "\n"
        )

        handle.write(
            f"Targets total: "
            f"{len(results)}\n"
        )

        handle.write(
            f"Targets PASS: "
            f"{pass_count}\n"
        )

        handle.write(
            f"Targets WARN: "
            f"{warn_count}\n"
        )

        handle.write(
            f"Targets failed: "
            f"{fail_count}\n"
        )

    print()
    print(
        "===== M03 SUMMARY ====="
    )

    print(
        f"Targets: "
        f"{len(results)}"
    )

    print(
        f"PASS: {pass_count}"
    )

    print(
        f"WARN: {warn_count}"
    )

    print(
        f"Failed: {fail_count}"
    )

    print(
        f"Summary: "
        f"{summary_path}"
    )


if __name__ == "__main__":
    main()
