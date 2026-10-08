#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import yaml


ROOT = Path.home() / "PocketToDynamics"

CONFIG = ROOT / "config/candidate_gate.yaml"

INPUT = (
    ROOT /
    "work/chemoinformatics/"
    "candidate_chemoinformatics.tsv"
)

OUT_DIR = ROOT / "work/candidates"

OUT_ALL = (
    OUT_DIR /
    "candidate_gate.tsv"
)

OUT_PROMOTED = (
    OUT_DIR /
    "promoted_candidates.tsv"
)

SUMMARY = (
    OUT_DIR /
    "candidate_gate_summary.json"
)

MODULE_LOG = (
    ROOT /
    "logs/M09_candidate_gate.log"
)


def load_tsv(path: Path):

    with path.open() as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


def as_int(value):

    return int(float(value))


def as_float(value):

    return float(value)


def main():

    if not INPUT.exists():

        print(
            "FAIL_M09_INPUT_NOT_FOUND"
        )

        sys.exit(1)

    cfg = yaml.safe_load(
        CONFIG.read_text()
    )

    gate = cfg["candidate_gate"]

    criteria = gate["criteria"]

    mode = (
        gate.get(
            "logic",
            {}
        ).get(
            "mode",
            "all",
        )
    )

    if mode not in {
        "all",
        "any",
    }:

        print(
            "FAIL_M09_INVALID_GATE_MODE"
        )

        sys.exit(1)

    rows = load_tsv(INPUT)

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    evaluated = []

    promoted = []

    for row in rows:

        checks = {
            "rank_pass":
                as_int(
                    row[
                        "integrated_rank"
                    ]
                )
                <= int(
                    criteria[
                        "max_integrated_rank"
                    ]
                ),

            "score_pass":
                as_float(
                    row[
                        "best_score"
                    ]
                )
                <= float(
                    criteria[
                        "max_best_score_kcal_mol"
                    ]
                ),

            "hit_fraction_pass":
                as_float(
                    row[
                        "hit_fraction"
                    ]
                )
                >= float(
                    criteria[
                        "min_hit_fraction"
                    ]
                ),

            "ligand_efficiency_pass":
                as_float(
                    row[
                        "docking_ligand_efficiency"
                    ]
                )
                >= float(
                    criteria[
                        "min_docking_ligand_efficiency"
                    ]
                ),
        }

        values = list(
            checks.values()
        )

        if mode == "all":

            gate_pass = all(
                values
            )

        else:

            gate_pass = any(
                values
            )

        reasons_failed = [
            key.replace(
                "_pass",
                ""
            )
            for key, value
            in checks.items()
            if not value
        ]

        evaluated_row = {
            "target_id":
                row["target_id"],

            "ligand_id":
                row["ligand_id"],

            "integrated_rank":
                row["integrated_rank"],

            "best_pocket":
                row["best_pocket"],

            "best_score":
                row["best_score"],

            "hit_fraction":
                row["hit_fraction"],

            "docking_ligand_efficiency":
                row[
                    "docking_ligand_efficiency"
                ],

            **{
                key:
                    int(value)
                for key, value
                in checks.items()
            },

            "gate_pass":
                int(gate_pass),

            "failed_criteria":
                ",".join(
                    reasons_failed
                ),
        }

        evaluated.append(
            evaluated_row
        )

        if gate_pass:

            promoted.append(
                evaluated_row
            )

        state = (
            "PROMOTE"
            if gate_pass
            else "REJECT"
        )

        print(
            f"[M09] "
            f"{row['target_id']} "
            f"{row['ligand_id']} "
            f"{state}"
        )

    fields = [
        "target_id",
        "ligand_id",
        "integrated_rank",
        "best_pocket",
        "best_score",
        "hit_fraction",
        "docking_ligand_efficiency",
        "rank_pass",
        "score_pass",
        "hit_fraction_pass",
        "ligand_efficiency_pass",
        "gate_pass",
        "failed_criteria",
    ]

    with OUT_ALL.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            evaluated
        )

    with OUT_PROMOTED.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            promoted
        )

    summary = {
        "evaluated":
            len(evaluated),

        "promoted":
            len(promoted),

        "rejected":
            len(evaluated)
            - len(promoted),

        "gate_mode":
            mode,

        "criteria":
            criteria,
    }

    SUMMARY.write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n"
    )

    if not evaluated:

        status = (
            "SKIP_M09_NO_CANDIDATES"
        )

    elif promoted:

        status = (
            "PASS_M09_CANDIDATE_GATE"
        )

    else:

        status = (
            "PASS_M09_NO_CANDIDATES_PROMOTED"
        )

    MODULE_LOG.write_text(
        "\n".join(
            [
                status,
                f"Evaluated: {len(evaluated)}",
                f"Promoted: {len(promoted)}",
                (
                    "Rejected: "
                    f"{len(evaluated)-len(promoted)}"
                ),
                f"Gate mode: {mode}",
                f"All candidates: {OUT_ALL}",
                (
                    "Promoted candidates: "
                    f"{OUT_PROMOTED}"
                ),
            ]
        )
        + "\n"
    )

    print()
    print(status)
    print(
        "Evaluated:",
        len(evaluated),
    )
    print(
        "Promoted:",
        len(promoted),
    )
    print(
        "Rejected:",
        len(evaluated)
        - len(promoted),
    )


if __name__ == "__main__":
    main()
