#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import hashlib
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

    v = as_float(value)
    if not v.is_integer() or v < 1:
        raise ValueError("Rank must be a positive integer")
    return int(v)


def as_float(value):

    v = float(value)
    if not math.isfinite(v):
        raise ValueError("Nonfinite candidate metric")
    return v


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

    as_int(criteria['max_integrated_rank'])
    for key, value in criteria.items():
        as_float(value)
    if not 0 <= as_float(criteria['min_hit_fraction']) <= 1:
        raise ValueError('Invalid hit-fraction threshold')
    if gate.get('output', {}).get('promote_failed_candidates', False):
        raise ValueError('Promoting failed candidates is unsupported')

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
    invalid_rows = 0

    for row in rows:

        invalid_data = False
        try:
            fraction = as_float(row['hit_fraction'])
            if not 0 <= fraction <= 1:
                raise ValueError('hit_fraction outside [0, 1]')
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
        except (ValueError, TypeError, KeyError, OverflowError):
            invalid_data = True
            invalid_rows += 1
            checks = dict.fromkeys(['rank_pass', 'score_pass', 'hit_fraction_pass', 'ligand_efficiency_pass'], False)

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

        gate_pass = gate_pass and not invalid_data

        reasons_failed = [
            key.replace(
                "_pass",
                ""
            )
            for key, value
            in checks.items()
            if not value
        ]

        if invalid_data:
            reasons_failed.append("invalid_input")

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
        "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "invalid_rows": invalid_rows,
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

    if invalid_rows:
        status = "FAIL_M09_INVALID_CANDIDATE_DATA"

    elif not evaluated:

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

    if invalid_rows:
        sys.exit(1)

if __name__ == "__main__":
    main()
