#!/usr/bin/env bash

set -euo pipefail

ROOT="$HOME/PocketToDynamics"

cd "$ROOT"

run_module() {

    local label="$1"
    shift

    echo
    echo "============================================================"
    echo "$label"
    echo "============================================================"

    "$@"
}


run_module \
    "M01 — Project initialization" \
    python scripts/01_init_project.py


run_module \
    "M02 — Input normalization" \
    python scripts/02_normalize_inputs.py


run_module \
    "M03 — Pocket mapping" \
    python scripts/03_map_pockets.py


run_module \
    "M04 — Job matrix" \
    python scripts/04_build_job_matrix.py


run_module \
    "M05A — Engine preflight" \
    python scripts/05_engine_preflight.py


run_module \
    "M05B — Molecular preparation" \
    python scripts/05_prepare_inputs.py


run_module \
    "M05C — Screening execution" \
    python scripts/05c_run_screening.py


run_module \
    "M06 — Screening parser" \
    python scripts/06_parse_screening.py


run_module \
    "M07 — Multi-pocket ranking" \
    python scripts/07_rank_candidates.py


run_module \
    "M08 — Chemoinformatics" \
    python scripts/08_chemoinformatics.py


run_module \
    "M09 — Candidate gate" \
    python scripts/09_candidate_gate.py


echo
echo "============================================================"
echo "PIPELINE STATUS"
echo "============================================================"

python scripts/update_status.py


echo
echo "============================================================"
echo "QC TESTS"
echo "============================================================"

python tests/test_m01_m04_integration.py
python tests/test_m06_parsing.py
python tests/test_m07_ranking.py
python tests/test_m08_chemoinformatics.py
python tests/test_m09_candidate_gate.py


echo
echo "============================================================"
echo "PASS_POCKET2DYNAMICS_CORE_PIPELINE"
echo "============================================================"
