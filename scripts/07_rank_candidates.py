#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import yaml


ROOT = Path.home() / "PocketToDynamics"

INPUT = ROOT / "work/parsed/task_summary.tsv"
CONFIG = ROOT / "config/ranking.yaml"

OUT_DIR = ROOT / "work/ranking"
OUT_TABLE = OUT_DIR / "candidate_ranking.tsv"

MODULE_LOG = ROOT / "logs/M07_ranking.log"


def load_tsv(path: Path):
    with path.open() as f:
        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


def mean(values):
    return statistics.mean(values) if values else math.nan


def median(values):
    return statistics.median(values) if values else math.nan


def stdev(values):
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def dense_rank(
    rows,
    field,
    reverse=False,
):

    values = sorted(
        {
            row[field]
            for row in rows
        },
        reverse=reverse,
    )

    lookup = {
        value: rank
        for rank, value in enumerate(
            values,
            start=1,
        )
    }

    for row in rows:
        row[
            f"rank_{field}"
        ] = lookup[row[field]]


def fmt(value):

    if value is None:
        return ""

    if isinstance(value, float):
        return f"{value:.6f}"

    return value


def main():

    if not INPUT.exists():

        print(
            "FAIL_M07_TASK_SUMMARY_NOT_FOUND"
        )

        sys.exit(1)

    cfg = yaml.safe_load(
        CONFIG.read_text()
    )

    hit_threshold = float(
        cfg["ranking"][
            "hit_threshold_kcal_mol"
        ]
    )

    weights = cfg["ranking"][
        "integrated_rank"
    ]

    task_rows = load_tsv(INPUT)

    grouped = defaultdict(list)

    for row in task_rows:

        key = (
            row["target_id"],
            row["ligand_id"],
        )

        grouped[key].append(row)

    candidates = []

    for (
        target_id,
        ligand_id,
    ), rows in sorted(grouped.items()):

        pocket_scores = [
            float(row["best_score"])
            for row in rows
        ]

        ordered = sorted(
            zip(
                pocket_scores,
                [
                    row["pocket_id"]
                    for row in rows
                ],
            ),
            key=lambda x: x[0],
        )

        best_score = ordered[0][0]
        best_pocket = ordered[0][1]

        if len(ordered) >= 2:

            pocket_selectivity_gap = (
                ordered[1][0]
                - ordered[0][0]
            )

        else:

            pocket_selectivity_gap = None

        n_hits = sum(
            score <= hit_threshold
            for score in pocket_scores
        )

        n_pockets = len(
            pocket_scores
        )

        hit_fraction = (
            n_hits / n_pockets
            if n_pockets
            else 0.0
        )

        candidates.append(
            {
                "target_id":
                    target_id,

                "ligand_id":
                    ligand_id,

                "n_pockets":
                    n_pockets,

                "n_hits":
                    n_hits,

                "hit_fraction":
                    hit_fraction,

                "best_score":
                    best_score,

                "mean_score":
                    mean(
                        pocket_scores
                    ),

                "median_score":
                    median(
                        pocket_scores
                    ),

                "score_sd":
                    stdev(
                        pocket_scores
                    ),

                "worst_score":
                    max(
                        pocket_scores
                    ),

                "best_pocket":
                    best_pocket,

                "pocket_selectivity_gap":
                    pocket_selectivity_gap,
            }
        )

    if not candidates:

        status = (
            "SKIP_M07_NO_CANDIDATES"
        )

        MODULE_LOG.write_text(
            status + "\n"
        )

        print(status)
        return

    # More negative = better
    dense_rank(
        candidates,
        "best_score",
        reverse=False,
    )

    dense_rank(
        candidates,
        "mean_score",
        reverse=False,
    )

    # Smaller SD = more robust
    dense_rank(
        candidates,
        "score_sd",
        reverse=False,
    )

    # Larger hit fraction = better
    dense_rank(
        candidates,
        "hit_fraction",
        reverse=True,
    )

    for row in candidates:

        row["integrated_rank_score"] = (

            row["rank_best_score"]
            * float(
                weights[
                    "rank_best_weight"
                ]
            )

            +

            row["rank_mean_score"]
            * float(
                weights[
                    "rank_mean_weight"
                ]
            )

            +

            row["rank_hit_fraction"]
            * float(
                weights[
                    "rank_hit_fraction_weight"
                ]
            )

            +

            row["rank_score_sd"]
            * float(
                weights[
                    "rank_robustness_weight"
                ]
            )
        )

    candidates.sort(
        key=lambda row: (
            row[
                "integrated_rank_score"
            ],
            row["best_score"],
            row["mean_score"],
            row["ligand_id"],
        )
    )

    # Final deterministic rank
    for rank, row in enumerate(
        candidates,
        start=1,
    ):

        row["integrated_rank"] = rank

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    fields = [
        "integrated_rank",
        "target_id",
        "ligand_id",
        "n_pockets",
        "n_hits",
        "hit_fraction",
        "best_score",
        "mean_score",
        "median_score",
        "score_sd",
        "worst_score",
        "best_pocket",
        "pocket_selectivity_gap",
        "rank_best_score",
        "rank_mean_score",
        "rank_hit_fraction",
        "rank_score_sd",
        "integrated_rank_score",
    ]

    with OUT_TABLE.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for row in candidates:

            writer.writerow(
                {
                    key: fmt(
                        row.get(key)
                    )
                    for key in fields
                }
            )

    status = (
        "PASS_M07_RANKING"
    )

    MODULE_LOG.write_text(
        "\n".join(
            [
                status,
                (
                    f"Candidates: "
                    f"{len(candidates)}"
                ),
                (
                    "Hit threshold: "
                    f"{hit_threshold}"
                ),
                (
                    f"Output: "
                    f"{OUT_TABLE}"
                ),
            ]
        )
        + "\n"
    )

    print()
    print(status)

    print(
        "Candidates:",
        len(candidates),
    )

    print(
        "Hit threshold:",
        hit_threshold,
    )

    print(
        "Output:",
        OUT_TABLE,
    )


if __name__ == "__main__":
    main()
